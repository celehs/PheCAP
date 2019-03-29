phecap_generate_dictionary_file <- function(
  cui_list, dict_file,
  user = "username", password = "password",
  host = "localhost", dbname = "umls", ...)
{
  if (length(cui_list) == 0L) {
    stop("'cui_list' is empty")
  }

  filter <- paste(paste0("'", cui_list, "'"), collapse = ", ")
  query <- sprintf("
    SELECT DISTINCT LOWER(str) AS str, cui
    FROM mrconso
    WHERE suppress='N' AND
    cui in (%s);
    ", filter)

  conn <- dbConnect(
    MySQL(),
    user = user,
    password = password,
    host = host,
    dbname = dbname,
    ...)

  rs <- dbSendQuery(conn, query)
  df <- dbFetch(rs)

  dbClearResult(rs)
  dbDisconnect(conn)

  if (nrow(df) == 0L) {
    stop("None of the CUIs are found in the database")
  }

  df <- df[!grepl("[;,\\[\\{\\(_<>=]", df$str, perl = TRUE), ]
  df <- df[!grepl(" - ", df$str, fixed = TRUE), ]
  df <- df[order(df$cui), ]

  write.table(
    df, dict_file,
    sep = "|", quote = FALSE,
    row.names = FALSE, col.names = FALSE)
  invisible(df)
}


phecap_perform_majority_voting <- function(
  input_folder)
{
  file_list <- list.files(input_folder, full.names = TRUE)
  if (length(file_list) == 0L) {
    stop("'input_folder' is empty or does not exist")
  }

  cui_each <- lapply(file_list, function(file) {
    textlines <- readLines(file)
    cui <- gregexpr("C[0-9]{7}", textlines, perl = TRUE)
    cui <- regmatches(textlines, cui)
    unique(unlist(cui))
  })

  cui <- table(unlist(cui_each))
  cui <- cui[cui >= length(file_list) * 0.5]
  cui <- sort.int(names(cui))

  cui
}
