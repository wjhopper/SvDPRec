packages <- c("dplyr","rprojroot")
loaded <- search()
for (p in packages) {
  if (! paste0("package:", p) %in% loaded) {
    library(p, character.only = TRUE)
  }
}
rm(packages, loaded, p)


load_data <- function(phase = c("Test","Slist","Points","Prac","TwoBack")) {
  
  phase <- tryCatch(match.arg(phase, several.ok = FALSE),
                    error = function(e) {
                      e$message <- sub(pattern = 'arg', "phase", e$message)
                      stop(e)
                    }
                    )
  root_dir <- rprojroot::is_rstudio_project$find_file()
  prefix <- paste0("RTBias_", phase, "_.+\\.txt$")
  file_names <- list.files(file.path(root_dir, "data"),
                           pattern = prefix,
                           full.names = TRUE)
  # dataset <- do.call(rbind,
  #                    lapply(file_names, read.table,
  #                           sep=",", stringsAsFactors=FALSE)
  #                    )
  dataset_list <- lapply(file_names, read.table,
                         sep=",", stringsAsFactors=FALSE)
  names(dataset_list) <- file_names
  dataset <- bind_rows(dataset_list, .id = "session")
  fname_length <- nchar(dataset$session)
  dataset$session <- as.numeric(substr(dataset$session, fname_length-4, fname_length-4))

  if (phase == "Test") {

    colnames(dataset) <- c("session","subject","list","trial", "pOld", "type", "strength",
                           "speeded_judgment", "speeded_RT", "V9",
                           "delayed_judgment", "delayed_RT", "V12",
                           "word", "V14")
    dataset <- dataset[setdiff(colnames(dataset), c("V9","V12","V14"))]
    dataset$type[dataset$word == "XXXXXXX"] <- NA_character_
    dataset$strength[dataset$strength == "bbb"] <- "L"
    dataset$subject <- as.character(dataset$subject)
  
  }
  
  return(dataset)
}