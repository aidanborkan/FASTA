# FASTA
FASTA Query for Proteomics

This script is designed for Rstudio and encapsulates a function to pull data from FASTA files from UniProt database.
A resulting dataframe will include select data based on specifications.

# Define a function that processes a FASTA file and finds shared protein entries
process_fasta_file <- function(fasta_file_path, comparison_df, output_file_name) {
  
  # Check for the comparison data frame and its 'id' column
  if (missing(comparison_df) || !"id" %in% names(comparison_df)) {
    stop("The comparison data frame is missing or does not have an 'id' column - rename applicable column.")
  }
  
  # Read the FASTA file
  fasta_data <- readAAStringSet(fasta_file_path)
  
  # Extract the headers
  headers <- names(fasta_data)
  
  # Extract UniProt names and gene names
  uniprot_names <- sub(".*\\|(.*?)\\|.*", "\\1", headers)
  collapsed_names <- sub(".*\\|.*?\\|(.*?)\\s.*", "\\1", headers)
  gene_names <- sub(".*GN=(\\w+).*", "\\1", headers)
  
  # Store results in a data.frame
  result_df <- data.frame(Header = headers, 
                          UniProt_Name = uniprot_names, 
                          Collapsed_Name = collapsed_names, 
                          Gene_Name = gene_names, 
                          stringsAsFactors = FALSE)
  
  # Find shared entries
  shared_entries <- intersect(result_df$UniProt_Name, comparison_df$id)
  
  # Print and write shared entries
  print(shared_entries)
  
  if (!missing(output_file_name)) {
    write.table(shared_entries, output_file_name, sep=",", row.names=FALSE, quote=FALSE)
  }
  
  # Return the result data frame and the shared entries
  return(list(Result_DataFrame = result_df, Shared_Entries = shared_entries))
}
# Invoke the function
#In this example the .gz FASTA file was downloaded from the UniProt database and unzipped
#In order: file path specification, modify Your_df with the name of the pre-existing dataframe, and the desired name for the result
shared_entries_cu <- process_fasta_file("uniprotkb_copper_binding_proteins_AND_m_2023_11_03.fasta.gz",Your_df, "Desired_Resulting_df.csv")
