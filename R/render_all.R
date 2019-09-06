library(rprojroot)

render_all <- function(report_dir="vignette", exclude="joint_diffusion.Rmd") {

  root_dir <- rprojroot::is_rstudio_project$find_file()
  
  Rmd_files <- list.files(file.path(root_dir, report_dir),
                          pattern = "\\.rmd$",
                          ignore.case = TRUE,
                          full.names = TRUE
                          )
  Rmd_files <- setdiff(Rmd_files, file.path(root_dir, report_dir, exclude))
  
  output_files <- rep(list(character(2)), length(Rmd_files))
  
  command_template <- 'Rscript -e "x <- rmarkdown::render(\\"%s\\", encoding=\\"UTF-8\\"); cat(x, sep="\n");"'
  
  for (i in 1:length(Rmd_files)) {

    f <- Rmd_files[i]
    message("Rendering file ", f)
    x <- system(sprintf(command_template, f), inter=TRUE, wait = TRUE)
    
    if (is.null(attr(x, "status"))) {
      len_x <- length(x)
      output_file <- x[len_x]
      cat(paste0(x[1:(len_x-1)], collapse = "\n"), sep="\n")
      output_files[[i]][1] <- output_file
    } else {
      cat(x)
      stop("Error rendering file", f)
    }
    
    files_dir <-  paste0(tools::file_path_sans_ext(output_file), "_files")
    if (file.exists(files_dir)) {
      output_files[[i]][2] <- files_dir
    }
  }
  
  return(output_files)
}


deploy_to_ghpages <- function(report_dir="vignette", exclude="joint_diffusion.Rmd") {
  
  root_dir <- rprojroot::is_rstudio_project$find_file()
  report_files <- render_all(report_dir, exclude)
  n_files_per_report <- vapply(report_files, FUN = length, FUN.VALUE = numeric(1))
  files_to_commit <- character(sum(n_files_per_report))
  
  insert_at <- 1
  for (i in 1:length(report_files)) {
    f <- report_files[[i]]
    file.copy(f[1], root_dir, overwrite = TRUE)
    files_to_commit[insert_at] <- shQuote(file.path(root_dir, basename(f[1])))
    insert_at <- insert_at + 1
    # check if report has accompnying _files/ directory
    if (f[2] != "") {
      file.copy(f[2], root_dir, overwrite = TRUE, recursive = TRUE)
      files_to_commit[insert_at] <- shQuote(file.path(root_dir, basename(f[2])))
      insert_at <- insert_at + 1
    }
  }

  sha <- system("git rev-parse --short HEAD", intern = TRUE, wait=TRUE)
  system("git checkout -b gh-pages")
  system(paste0("git add ", files_to_commit, collapse = " "), wait = TRUE)
  system(sprintf("git commit -m 'Updating from %s'", sha), wait = TRUE)
    
}

