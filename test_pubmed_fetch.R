# Smoke test for the new PubMed fetch logic in cv.qmd.
# Run interactively:  source("test_pubmed_fetch.R")
#
# This pulls the same fetch_myncbi() / fetch_rentrez() pipeline the CV uses,
# then prints the first 5 numbered citations so you can eyeball the format
# before re-rendering cv.qmd.

suppressPackageStartupMessages({
  library(yaml)
  library(rvest)
  library(xml2)
  library(rentrez)
})

cv_data    <- read_yaml("cv_data.yml")
myncbi_url <- cv_data$personal$links$pubmed

# ---- fetch_myncbi (paginated) ----
fetch_myncbi <- function(url, max_pages = 20) {
  safe_text <- function(node) {
    if (inherits(node, "xml_missing")) return("")
    txt <- rvest::html_text2(node)
    if (is.null(txt) || is.na(txt)) "" else txt
  }

  parse_page <- function(page) {
    items <- page %>% rvest::html_elements("div.citation")
    if (length(items) == 0) return(NULL)

    rows <- lapply(items, function(item) {
      docsum <- rvest::html_element(item, "div.ncbi-docsum")
      if (inherits(docsum, "xml_missing")) docsum <- item

      authors <- sub("\\.\\s*$", "", trimws(safe_text(
        rvest::html_element(docsum, ".authors"))))
      title_node <- rvest::html_element(docsum, "a[href*='/pubmed/']")
      title   <- trimws(safe_text(title_node))
      journal <- trimws(safe_text(rvest::html_element(docsum, ".source")))
      pubdate <- sub(";\\s*$", "", trimws(safe_text(
        rvest::html_element(docsum, ".pubdate"))))
      volume  <- trimws(safe_text(rvest::html_element(docsum, ".volume")))
      issue   <- trimws(safe_text(rvest::html_element(docsum, ".issue")))

      date <- pubdate
      if (nzchar(volume)) date <- paste0(date, ";", volume)
      if (nzchar(issue))  date <- paste0(date, issue)
      date <- trimws(date)

      pmid <- ""
      check <- rvest::html_element(item, "input.citation-check")
      if (!inherits(check, "xml_missing")) {
        v <- rvest::html_attr(check, "pmid")
        if (!is.null(v) && !is.na(v)) pmid <- v
      }
      if (!nzchar(pmid)) {
        ptxt <- safe_text(rvest::html_element(docsum, ".pmid"))
        m <- regmatches(ptxt, regexpr("\\d+", ptxt))
        if (length(m)) pmid <- m
      }
      if (!nzchar(pmid) && !inherits(title_node, "xml_missing")) {
        href <- rvest::html_attr(title_node, "href")
        if (!is.null(href) && !is.na(href)) {
          m <- regmatches(href, regexpr("\\d+", href))
          if (length(m)) pmid <- m
        }
      }

      data.frame(authors = authors, title = title, journal = journal,
                 date = date, pmid = pmid, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }

  tryCatch({
    # Strip any existing query string — MyNCBI ignores ?page=N when appended
    # to ?sort=date&direction=descending and just returns page 1.
    base_url <- sub("\\?.*$", "", url)
    acc <- list(); seen_pmids <- character(0)

    for (p in seq_len(max_pages)) {
      page_url <- paste0(base_url, "?page=", p)
      page <- rvest::read_html(page_url)
      df   <- parse_page(page)
      cat(sprintf("  page %d: %s citations\n", p,
                  if (is.null(df)) "0" else as.character(nrow(df))))
      if (is.null(df) || nrow(df) == 0) break

      new_rows <- df[!df$pmid %in% seen_pmids, , drop = FALSE]
      if (nrow(new_rows) == 0) break
      acc[[length(acc) + 1]] <- new_rows
      seen_pmids <- c(seen_pmids, new_rows$pmid)

      if (nrow(df) < 50) break
    }

    if (length(acc) == 0) return(NULL)
    df <- do.call(rbind, acc)
    df <- df[nchar(df$title) > 0 & nchar(df$pmid) > 0, , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df
  }, error = function(e) {
    message("MyNCBI scrape failed: ", conditionMessage(e)); NULL
  })
}

# ---- fetch_rentrez ----
fetch_rentrez <- function(orcid = "0000-0002-3341-9022",
                          author_term = "Hawken S[Author]",
                          aff_terms = c("Ottawa", "OHRI",
                                        "Institute for Clinical Evaluative Sciences",
                                        "ICES",
                                        "CHEO",
                                        "Children's Hospital of Eastern Ontario",
                                        "McMaster University",
                                        "Hamilton Health Sciences",
                                        "Population Health Research Institute"),
                          retmax = 1000) {
  tryCatch({
    aff <- paste0("(", paste(sprintf("%s[Affiliation]", aff_terms), collapse = " OR "), ")")
    aff_clause <- paste0("(", author_term, " AND ", aff, ")")
    query <- if (!is.null(orcid) && nzchar(orcid)) {
      paste0("(", orcid, "[AUID] OR ", aff_clause, ")")
    } else aff_clause
    cat("rentrez query: ", query, "\n", sep = "")
    sres <- rentrez::entrez_search(db = "pubmed", term = query, retmax = retmax)
    cat("rentrez found ", length(sres$ids), " PMIDs\n", sep = "")
    if (length(sres$ids) == 0) return(NULL)

    id_chunks <- split(sres$ids, ceiling(seq_along(sres$ids) / 100))
    xml_docs <- lapply(id_chunks, function(ids) {
      raw <- rentrez::entrez_fetch(db = "pubmed", id = ids,
                                   rettype = "xml", parsed = FALSE)
      xml2::read_xml(raw)
    })

    parse_article <- function(art) {
      gx <- function(xp) {
        n <- xml2::xml_find_first(art, xp)
        if (inherits(n, "xml_missing")) "" else xml2::xml_text(n)
      }
      pmid <- gx(".//PMID"); title <- gx(".//ArticleTitle")
      journal <- gx(".//Journal/Title")
      year <- gx(".//JournalIssue/PubDate/Year")
      if (!nzchar(year)) year <- substr(gx(".//JournalIssue/PubDate/MedlineDate"), 1, 4)
      volume <- gx(".//JournalIssue/Volume")
      issue  <- gx(".//JournalIssue/Issue")
      pages  <- gx(".//Pagination/MedlinePgn")

      authors_nodes <- xml2::xml_find_all(art, ".//AuthorList/Author")
      authors <- vapply(authors_nodes, function(a) {
        last <- xml2::xml_text(xml2::xml_find_first(a, "LastName"))
        init <- xml2::xml_text(xml2::xml_find_first(a, "Initials"))
        if (is.na(last)) {
          coll <- xml2::xml_text(xml2::xml_find_first(a, "CollectiveName"))
          if (is.na(coll)) "" else coll
        } else paste0(last, if (!is.na(init) && nzchar(init)) paste0(" ", init))
      }, character(1))
      authors <- paste(authors[nzchar(authors)], collapse = ", ")

      date <- year; vip <- ""
      if (nzchar(volume)) vip <- volume
      if (nzchar(issue))  vip <- paste0(vip, "(", issue, ")")
      if (nzchar(pages))  vip <- paste0(vip, ":", pages)
      if (nzchar(vip))    date <- paste0(date, ";", vip)

      data.frame(authors = authors, title = title, journal = journal,
                 date = date, pmid = pmid,
                 year = suppressWarnings(as.integer(substr(year, 1, 4))),
                 stringsAsFactors = FALSE)
    }

    all_rows <- list()
    for (doc in xml_docs) {
      arts <- xml2::xml_find_all(doc, ".//PubmedArticle")
      for (a in arts) all_rows[[length(all_rows) + 1]] <- parse_article(a)
    }
    df <- do.call(rbind, all_rows)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df <- df[order(-df$year, df$title), ]; df$year <- NULL; df
  }, error = function(e) {
    message("rentrez fetch failed: ", conditionMessage(e)); NULL
  })
}

# ---- Pagination probe: see what URLs MyNCBI uses for paging ----
cat("Probing MyNCBI page 1 for pagination links...\n")
tryCatch({
  p1 <- rvest::read_html(myncbi_url)
  links <- p1 %>% rvest::html_elements("a")
  hrefs <- rvest::html_attr(links, "href")
  paging <- unique(hrefs[grepl("page|start|offset|currentPage|p=", hrefs,
                               ignore.case = TRUE)])
  paging <- paging[!is.na(paging)]
  if (length(paging) > 0) {
    cat("Candidate paging hrefs found on page 1:\n")
    for (h in head(paging, 15)) cat("  ", h, "\n")
  } else {
    cat("No obvious pagination hrefs found in <a> tags.\n")
  }

  # Also report the citation count badge if present
  count_node <- p1 %>% rvest::html_element(".result-count, .citations-count, .total")
  if (!inherits(count_node, "xml_missing")) {
    cat("Count badge text: ", rvest::html_text2(count_node), "\n")
  }
}, error = function(e) cat("Probe failed: ", conditionMessage(e), "\n"))
cat("\n")

cat("Querying PubMed via rentrez (MyNCBI scrape disabled — JS pagination)...\n")
pubs <- fetch_rentrez()
if (is.null(pubs)) {
  cat("rentrez returned nothing — running MyNCBI diagnostic for the record...\n\n")
  tryCatch({
    page <- rvest::read_html(myncbi_url)
    cat("Page title: ", rvest::html_text2(rvest::html_element(page, "title")), "\n")
    html_len <- nchar(as.character(page))
    cat("HTML length: ", html_len, " chars\n")

    # Count candidate containers
    sels <- c("div.ncbi-docsum", "li.ncbi-docsum", "[data-pmid]",
              "li.citation", "div.citation", "article", "li.item",
              ".bibliography li", ".citation-list li")
    cat("\nSelector matches:\n")
    for (s in sels) {
      n <- length(rvest::html_elements(page, s))
      cat(sprintf("  %-30s  %d\n", s, n))
    }

    # Any data-pmid attributes anywhere?
    pmid_nodes <- rvest::html_elements(page, "[data-pmid]")
    if (length(pmid_nodes) > 0) {
      cat("\nFound", length(pmid_nodes), "nodes with data-pmid; first tag:\n")
      print(pmid_nodes[[1]])
    }

    # Dump a trimmed copy of the page for inspection
    dump_path <- "_myncbi_page_dump.html"
    writeLines(as.character(page), dump_path)
    cat("\nFull page written to ", dump_path,
        " — open it and search for one of your paper titles,\n",
        "then look at the surrounding <div>/<li> tags and class names.\n",
        "Send those class names back and I'll fix the selectors.\n", sep = "")
  }, error = function(e) {
    cat("Diagnostic failed too: ", conditionMessage(e), "\n")
  })
  cat("\nFalling back to MyNCBI scrape (will likely return only first 50)...\n")
  pubs <- fetch_myncbi(myncbi_url)
}

if (is.null(pubs) || nrow(pubs) == 0) {
  stop("Could not retrieve publications from MyNCBI or PubMed.")
}

cat(sprintf("\nRetrieved %d publications. First 5 as they'd appear in the CV:\n\n",
            nrow(pubs)))

bold_author <- function(text) {
  gsub("(Hawken S|Hawken, S|Steven Hawken|S Hawken|S\\.? Hawken)",
       "**\\1**", text, ignore.case = FALSE)
}

for (i in seq_len(min(5, nrow(pubs)))) {
  p <- pubs[i, ]
  line <- paste0(i, ". ", bold_author(p$authors),
                 ". **", sub("\\.\\s*$", "", p$title), "**. ",
                 if (nzchar(p$journal)) paste0("*", p$journal, "*. ") else "",
                 if (nzchar(p$date)) paste0(sub("\\.\\s*$", "", p$date), ". ") else "",
                 if (nzchar(p$pmid)) sprintf("[PMID: %s]", p$pmid) else "")
  cat(line, "\n\n")
}
