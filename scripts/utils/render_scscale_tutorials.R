#!/usr/bin/env Rscript

docs_dir <- "packages/scScale/docs"

render_one <- function(input) {
  old <- setwd(docs_dir)
  on.exit(setwd(old), add = TRUE)

  md <- knitr::knit(input, output = sub("[.]Rmd$", ".md", input), quiet = TRUE)
  raw <- readLines(input, warn = FALSE)
  title_line <- grep("^title:", raw, value = TRUE)[1]
  title <- sub("^title:[[:space:]]*", "", title_line)
  title <- gsub('^"|"$', "", title)

  lines <- readLines(md, warn = FALSE)
  if (length(lines) >= 3 && identical(lines[1], "---")) {
    ends <- which(lines[-1] == "---")
    if (length(ends)) lines <- lines[-seq_len(ends[1] + 1)]
  }

  body <- commonmark::markdown_html(paste(lines, collapse = "\n"), extensions = TRUE)
  html <- paste0(
    '<!doctype html>\n',
    '<html lang="en">\n',
    '<head>\n',
    '<meta charset="utf-8">\n',
    '<meta name="viewport" content="width=device-width, initial-scale=1">\n',
    '<title>', title, '</title>\n',
    '<style>',
    'body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;line-height:1.55;max-width:920px;margin:32px auto;padding:0 20px;color:#202124}',
    'pre{background:#f6f8fa;padding:12px;overflow:auto;border-radius:6px}',
    'code{font-family:Menlo,Consolas,monospace}',
    'img{max-width:100%;height:auto}',
    'table{border-collapse:collapse}',
    'td,th{border:1px solid #ddd;padding:4px 8px}',
    'blockquote{color:#555;border-left:4px solid #ddd;margin-left:0;padding-left:12px}',
    'h1{margin-bottom:0.4em}',
    '</style>\n',
    '</head>\n',
    '<body>\n',
    '<h1>', title, '</h1>\n',
    body,
    '\n</body>\n',
    '</html>\n'
  )
  writeLines(html, sub("[.]Rmd$", ".html", input))
}

files <- c(
  "scScale-tutorial.Rmd",
  "cell-number-scaling-law.Rmd",
  "umi-scaling-law.Rmd",
  "batch-number-scaling-law.Rmd"
)

invisible(lapply(files, render_one))
cat("Rendered", length(files), "scScale tutorials.\n")
