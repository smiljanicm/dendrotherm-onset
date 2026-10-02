library(zip)

folder_path <- "./output"
zip_file <- "./output/tables.zip"

all_items <- list.files(folder_path, full.names = TRUE)
files_only <- all_items[!file.info(all_items)$isdir]

zip::zip(
  zipfile = zip_file,
  files = files_only
)

folder_path <- "./output/german"
zip_file <- "./output/tables_german.zip"

zip::zip(
  zipfile = zip_file,
  files = list.files(folder_path, full.names = TRUE)
)