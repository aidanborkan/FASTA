# ------------------------------------------------------------
# FASTA_R Shiny
# UniProt FASTA Query + Header Parsing + Optional Overlap
#
# What it does:
# 1) User enters a UniProt query (e.g., "proteome:UP000005640 AND reviewed:true")
# 2) App fetches FASTA via UniProt REST stream endpoint
# 3) Parses FASTA headers into a tidy table: accession / entry name / gene
# 4) Optional: user uploads a CSV/TSV with protein IDs and computes overlap
#
# UniProt REST stream endpoint pattern:
# https://rest.uniprot.org/uniprotkb/stream?query=...&format=fasta&compressed=...
# (See UniProt programmatic access + API query docs.) :contentReference[oaicite:2]{index=2}
# ------------------------------------------------------------

options(shiny.maxRequestSize = 1000 * 1024^2)  # allow big uploads (~1GB)

library(shiny)
library(httr2)
library(Biostrings)
library(DT)

# ----------------------------
# Helpers
# ----------------------------

# Parse a UniProt FASTA header like:
#   sp|P12345|ENTRY_NAME Some protein OS=... GN=GENE ...
#   tr|Q9XYZ1|ENTRY_NAME ...
parse_uniprot_header <- function(header) {
  
  # Accessions are typically the token between first and second "|"
  accession <- if (grepl("\\|", header)) sub("^.*\\|([^|]+)\\|.*$", "\\1", header) else NA_character_
  
  # Entry name typically token after second "|", up to first space
  entry_name <- if (grepl("\\|", header)) sub("^.*\\|[^|]+\\|([^ ]+).*$", "\\1", header) else NA_character_
  
  # Gene name (GN=) is optional
  gene_name <- if (grepl("GN=", header)) sub("^.*GN=([^ ]+).*$", "\\1", header) else NA_character_
  
  data.frame(
    Header = header,
    UniProt_Accession = accession,
    UniProt_EntryName = entry_name,
    Gene_Name = gene_name,
    stringsAsFactors = FALSE
  )
}

# Download FASTA from UniProt stream endpoint to a tempfile and return path
fetch_uniprot_fasta <- function(query, include_isoform = FALSE, compressed = FALSE, user_agent = "FASTA_R_Shiny/1.0") {
  
  base <- "https://rest.uniprot.org/uniprotkb/stream"
  
  # UniProt supports query + format=fasta + compressed=... and includeIsoform=true/false
  # (Parameter names and endpoint pattern documented by UniProt and widely used.) :contentReference[oaicite:3]{index=3}
  params <- list(
    query = query,
    format = "fasta",
    compressed = if (compressed) "true" else "false",
    includeIsoform = if (include_isoform) "true" else "false"
  )
  
  # File target
  tmp <- tempfile(fileext = if (compressed) ".fasta.gz" else ".fasta")
  
  req <- request(base) |>
    req_url_query(!!!params) |>
    req_user_agent(user_agent) |>
    req_error(is_error = function(resp) resp_status(resp) >= 400)
  
  # Stream to disk
  resp <- req_perform(req)
  writeBin(resp_body_raw(resp), tmp)
  tmp
}

# Read FASTA file into AAStringSet and parse headers to a data.frame
parse_fasta_file <- function(fasta_path) {
  aa <- readAAStringSet(fasta_path)
  headers <- names(aa)
  result_df <- do.call(rbind, lapply(headers, parse_uniprot_header))
  list(aa = aa, parsed = result_df)
}

# Read comparison file (CSV/TSV) and return data.frame
read_comparison <- function(path) {
  # Try TSV first if it looks tabby; otherwise CSV
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

# ----------------------------
# UI
# ----------------------------
ui <- fluidPage(
  titlePanel("FASTA_R — UniProt FASTA Query (Proteomics)"),
  
  sidebarLayout(
    sidebarPanel(
      tags$h4("1) UniProt Query"),
      
      textInput(
        "uniprot_query",
        "UniProt query",
        value = "proteome:UP000005640 AND reviewed:true"
      ),
      helpText("Tip: UniProt query syntax supports fields like proteome:, organism_id:, reviewed:true, etc."),
      
      checkboxInput("include_isoform", "Include isoforms", value = FALSE),
      checkboxInput("compressed", "Download compressed (.gz)", value = FALSE),
      
      actionButton("fetch_btn", "Fetch FASTA from UniProt"),
      
      tags$hr(),
      
      tags$h4("2) Optional overlap vs your IDs"),
      fileInput("comparison_file", "Upload comparison IDs (CSV/TSV)", accept = c(".csv", ".tsv", ".txt")),
      uiOutput("id_column_ui"),
      
      selectInput(
        "match_on",
        "Match on FASTA identifier",
        choices = c("UniProt accession" = "UniProt_Accession",
                    "UniProt entry name" = "UniProt_EntryName",
                    "Gene name" = "Gene_Name"),
        selected = "UniProt_Accession"
      ),
      
      tags$hr(),
      
      downloadButton("download_fasta", "Download FASTA"),
      downloadButton("download_overlap", "Download overlap IDs")
    ),
    
    mainPanel(
      tags$h4("Status"),
      verbatimTextOutput("status"),
      
      tags$h4("Parsed FASTA headers"),
      DTOutput("parsed_table"),
      
      tags$h4("Overlap (if comparison file provided)"),
      DTOutput("overlap_table")
    )
  )
)

# ----------------------------
# Server
# ----------------------------
server <- function(input, output, session) {
  
  rv <- reactiveValues(
    fasta_path = NULL,
    parsed_df = NULL,
    aa = NULL,
    comparison_df = NULL,
    overlap_ids = character(0)
  )
  
  output$status <- renderText("Ready.")
  
  # Load comparison file if provided
  observeEvent(input$comparison_file, {
    req(input$comparison_file)
    df <- tryCatch(read_comparison(input$comparison_file$datapath),
                   error = function(e) NULL)
    if (is.null(df) || !is.data.frame(df) || ncol(df) < 1) {
      rv$comparison_df <- NULL
      output$status <- renderText("Comparison file could not be read (expect CSV/TSV).")
    } else {
      rv$comparison_df <- df
      output$status <- renderText(sprintf("Loaded comparison file with %d rows, %d columns.", nrow(df), ncol(df)))
    }
  })
  
  # Dynamic UI for choosing the ID column from comparison file
  output$id_column_ui <- renderUI({
    if (is.null(rv$comparison_df)) return(NULL)
    cols <- names(rv$comparison_df)
    selectInput("id_column", "ID column in comparison file", choices = cols, selected = cols[1])
  })
  
  # Fetch + parse FASTA
  observeEvent(input$fetch_btn, {
    req(input$uniprot_query)
    q <- trimws(input$uniprot_query)
    if (!nzchar(q)) {
      output$status <- renderText("Query is empty.")
      return()
    }
    
    output$status <- renderText("Fetching FASTA from UniProt...")
    
    fasta_path <- tryCatch(
      fetch_uniprot_fasta(
        query = q,
        include_isoform = isTRUE(input$include_isoform),
        compressed = isTRUE(input$compressed)
      ),
      error = function(e) {
        output$status <- renderText(paste("Fetch failed:", e$message))
        return(NULL)
      }
    )
    if (is.null(fasta_path)) return()
    
    parsed <- tryCatch(
      parse_fasta_file(fasta_path),
      error = function(e) {
        output$status <- renderText(paste("Parse failed:", e$message))
        return(NULL)
      }
    )
    if (is.null(parsed)) return()
    
    rv$fasta_path <- fasta_path
    rv$aa <- parsed$aa
    rv$parsed_df <- parsed$parsed
    
    output$status <- renderText(sprintf("Fetched and parsed FASTA. Entries: %d", nrow(rv$parsed_df)))
  })
  
  # Show parsed headers
  output$parsed_table <- renderDT({
    req(rv$parsed_df)
    datatable(rv$parsed_df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Compute overlap whenever relevant inputs change
  observe({
    req(rv$parsed_df)
    
    # Need comparison df + id column to compute overlap
    if (is.null(rv$comparison_df) || is.null(input$id_column) || !nzchar(input$id_column)) {
      rv$overlap_ids <- character(0)
      return()
    }
    
    fasta_ids <- rv$parsed_df[[input$match_on]]
    comp_ids  <- rv$comparison_df[[input$id_column]]
    
    fasta_ids <- trimws(as.character(fasta_ids))
    comp_ids  <- trimws(as.character(comp_ids))
    
    fasta_ids <- fasta_ids[!is.na(fasta_ids) & nzchar(fasta_ids)]
    comp_ids  <- comp_ids[!is.na(comp_ids) & nzchar(comp_ids)]
    
    rv$overlap_ids <- sort(unique(intersect(fasta_ids, comp_ids)))
  })
  
  output$overlap_table <- renderDT({
    if (length(rv$overlap_ids) == 0) {
      return(datatable(data.frame(Overlap_IDs = character(0)), options = list(pageLength = 5)))
    }
    datatable(data.frame(Overlap_IDs = rv$overlap_ids), options = list(pageLength = 10))
  })
  
  # Download FASTA
  output$download_fasta <- downloadHandler(
    filename = function() {
      if (isTRUE(input$compressed)) "uniprot_query.fasta.gz" else "uniprot_query.fasta"
    },
    content = function(file) {
      req(rv$fasta_path)
      file.copy(rv$fasta_path, file, overwrite = TRUE)
    }
  )
  
  # Download overlap IDs
  output$download_overlap <- downloadHandler(
    filename = function() "overlap_ids.csv",
    content = function(file) {
      write.csv(data.frame(Overlap_IDs = rv$overlap_ids), file, row.names = FALSE)
    }
  )
  
}

shinyApp(ui = ui, server = server)
