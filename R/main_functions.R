#' @docType package
#' @name bedanno
#' @title Annotate BED-like variant files.
#' @import dplyr
NULL
#' @import data.table
NULL
#' @import doMC
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to bedanno! Both variants and annotations files should follow the BED schema for the first three columns.")
}

# Utility Functions -------------------------------------------------

#' Read gzipped file into data.table.
#'
#' @family utility functions
#' @title gzfread
#' @description Import text from gzipped or gunzipped file to data.table
#' @aliases gzfread
#' @author http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' @param path a character
#' @return data.table
gzfread <- function(path, sep = "\t", out_dir = NULL){
  if(!stringr::str_detect(path, ".gz$")) {
    return(data.table::fread(path, sep = sep))
  }
  else{
    return(data.table::fread(paste0("zcat < ", "\"", path, "\"")))
  }
}

#' Import BED file.
#'
#' @family utility functions
#' @param bed_path A character.
#' @return data.table
get_bed <- function(bed_path){
  paste("Reading ", basename(bed_path), sep="") %>% message
  this_bed <- gzfread(bed_path, sep="\t")
  # Alternative:
  # this_bed <- read_delim(bed_path, delim="\t", col_names=FALSE)
  data.table::setnames(this_bed, names(this_bed)[1:3],
                       c("CHROM", "START", "STOP")) %>%
    select(CHROM, START, STOP) %>%
    as.data.table %>%
    return
}

#' Import variant table in pseudo-BED format.
#'
#' @family utility functions
#' @param variant_path A character.
#' @return data.table
get_variant_table <- function(variant_path){
  variants <- gzfread(variant_path, sep="\t")
  data.table::setnames(variants, names(variants)[1:3], c("CHROM", "START", "STOP")) %>%
    data.table::setkeyv(c("CHROM", "START", "STOP")) %>%
    # Vitally important to get unique rows by CHROM, START, STOP; way to keep
    # annotation rows aligned with data table rows (must always keep same sort)
    unique %>%
    return
}

#' Get the full paths of files from a directory.
#'
#' @family utility functions
#' @param dir_path A character.
#' @return character
get_file_paths <- function(dir_path){
  file_list <- paste(dir_path, "/", list.files(dir_path, recursive = TRUE), sep="")
  return(file_list)
}

# Main Functions ----------------------------------------------------

#' Generate binary annotation column.
#'
#' @family annotation functions with NO intermediate file output
#' @param bed_file_path A character.
#' @param variant_table A data.table.
#' @return Annotation column.
get_anno_col <- function(bed_file_path, variant_table){
  bed <- get_bed(bed_file_path)
  try_output <- try({
    data.table::setkey(bed, CHROM, START, STOP)
    overlap_table <- foverlaps(variant_table, bed, type="within")
    distinct_overlap_table <- unique(overlap_table, by=c("CHROM", "i.START", "i.STOP"))
    anno_column <- rep.int(0, nrow(distinct_overlap_table)) # Use rep.int instead of rep for speed
    anno_column[!is.na(distinct_overlap_table$START)] <- 1
    anno_column
  }, silent = TRUE)

  # If there is something wrong with the BED file return a column of NAs
  if(class(try_output) == "try-error"){
    print(paste0("BED annotation file ", bed_file_path, " is malformed. Inserting NAs."))
    return(rep(NA, nrow(unique(variant_table, by = c("CHROM", "START", "STOP")))))
  }else{
    return(anno_column)
  }
}

#' Annotate a data.table with binary annotations.
#'
#' @family annotation functions with NO intermediate file output
#' @param bed_dir_path A character.
#' @param variant_path A data.table.
#' @param cores A numeric.
#' @return Annotated data.frame.
annotate_variants <- function(bed_dir_path, variant_path, cores=1){
  doMC::registerDoMC(cores=cores)
  elapsed_get_variant_table <- system.time(variant_table <- get_variant_table(variant_path))
  message(paste("Variant table read in", elapsed_get_variant_table["elapsed"], "seconds."))
  bed_paths <- get_file_paths(bed_dir_path)
  # Old Method:
  # elapsed_anno_list <- (system.time(anno_list <- parallel::mclapply(bed_paths, get_anno_col, variant_table, mc.cores=cores)))
  elapsed_anno_list <- system.time(anno_list <- plyr::alply(bed_paths, 1, get_anno_col, variant_table, .parallel = TRUE))
  message(paste("Annotation matrix generated in", elapsed_anno_list["elapsed"], "seconds"))
  names(anno_list) <- basename(bed_paths)
  doMC::registerDoMC(cores=1)

  anno_df <- data.frame(variant_table, anno_list)
  write.table(anno_df, file = paste0(variant_path, ".anno", sep=""), col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")
  return(anno_df)
}
# -------------------------------------------------------------------

#' Write intermediate annotation column.
#'
#' @family annotation functions WITH intermediate file output
#' @param bed_file_path A character.
#' @param variant_table A data.table.
#' @param anno_path A character.
#' @return NULL.
write_anno_col <- function(bed_file_path, variant_table, anno_path){
  bed <- get_bed(bed_file_path)
  try_output <- try({
    data.table::setkey(bed, CHROM, START, STOP)
    overlap_table <- foverlaps(variant_table, bed, type="within")
    distinct_overlap_table <- unique(overlap_table, by=c("CHROM", "i.START", "i.STOP"))
    anno_column <- rep.int(0, nrow(distinct_overlap_table)) # Use rep.int instead of rep for speed
    anno_column[!is.na(distinct_overlap_table$START)] <- 1
    anno_column <- data.frame(anno_column)
    names(anno_column) <- basename(bed_file_path)
    }, silent = TRUE)

  # If there is something wrong with the BED file print a column of NAs
  if(class(try_output) == "try-error"){
    print(paste0("BED annotation file ", bed_file_path, " is malformed. Inserting NAs."))
    na_column <- rep(NA, nrow(unique(variant_table, by = c("CHROM", "START", "STOP"))))
    na_column <- data.frame(na_column)
    names(na_column) <- basename(bed_file_path)
    write.table(na_column, file = paste(anno_path, "/", basename(bed_file_path), ".anno", sep=""), col.names = TRUE, quote = FALSE, row.names = FALSE)}else{
    write.table(anno_column, file = paste(anno_path, "/", basename(bed_file_path), ".anno", sep=""), col.names = TRUE, quote = FALSE, row.names = FALSE)
    }
}

#' Write intermediate annotation columns from dir.
#'
#' @family annotation functions WITH intermediate file output
#' @param bed_dir_path A character.
#' @param variant_path A data.table.
#' @param cores A numeric.
#' @return NULL.
write_variant_annotations <- function(bed_dir_path, variant_path, cores = 1){
  anno_dir_path <- file.path(dirname(variant_path), paste(basename(variant_path), "annos", sep="."))
  dir.create(anno_dir_path, showWarnings = FALSE)
  elapsed_get_variant_table <- system.time(variant_table <- get_variant_table(variant_path))
  write.table(variant_table, paste(variant_path, ".reordered", sep=""), col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")
  message(paste("Variant table read in", elapsed_get_variant_table["elapsed"], "seconds."))
  bed_paths <- get_file_paths(bed_dir_path)
  elapsed_anno <- (system.time(parallel::mclapply(bed_paths, write_anno_col, variant_table, anno_dir_path, mc.cores = cores)))
  message(paste("Annotation files written in", elapsed_anno["elapsed"], "seconds"))
  message("Cleaning up intermediate gunzipped BED files.")
  new_bed_paths <- get_file_paths(bed_dir_path)
  file.remove(new_bed_paths[!(new_bed_paths %in% bed_paths)])
}

#' Horizontally concatenate intermediate annotation columns.
#'
#' @family annotation functions WITH intermediate file output
#' @param variant_path A character.
#' @param anno_col_dir_path A character.
#' @return NULL.
horizontal_concat_annos <- function(variant_path, anno_col_dir_path){
  anno_file_names <- (list.files(anno_col_dir_path)) %>% stringr::str_detect(".anno$") %>% list.files(anno_col_dir_path)[.]
  paths <- paste(anno_col_dir_path, anno_file_names, sep="/")
  paths <- paste(paths, collapse = " ")
  paths <- paste(variant_path, paths, collapse = " ")
  # Use UNIX paste command; set soft limit to open number of files to maximum (10240)
  sys_command <- paste("ulimit -Sn 10240; paste ", paths, " > ", variant_path, ".annotated", sep="")
  # message("Horizontal concatentation system command:")
  message(sys_command)
  message(paste("Adding", length(unlist(stringr::str_split(paths, " "))) - 1, "annotation columns to variant file."))
  system(sys_command)
}

#' Write intermediate cols and horizontally concatenate output.
#'
#' @family annotation functions WITH intermediate file output
#' @param bed_dir_path A character.
#' @param variant_path A data.table.
#' @param cores A numeric.
#' @return NULL.
annotate_variants_with_intermediates <- function(bed_dir_path, variant_path, cores = 1){
  write_variant_annotations(bed_dir_path = bed_dir_path, variant_path = variant_path, cores = cores)
  anno_col_dir_path <- file.path(dirname(variant_path), paste(basename(variant_path), "annos", sep="."))
  formatted_var_path <- paste(variant_path, ".reordered", sep="")
  horizontal_concat_annos(variant_path = formatted_var_path, anno_col_dir_path = anno_col_dir_path)
}

# Top-level functions -----------------------------------------------
#' Annotate variant file from BED files
#'
#' @family Annotate variant file with returned data.table or annotated files.
#' @param bed_dir_path A character.
#' @param variant_path A data.table.
#' @param cores A numeric.
#' @param write_files A logical.
#' @return NULL for write_files == FALSE; data.table for write_files == TRUE
#' @export
anno_vars <- function(bed_dir_path, variant_path, cores = 1, write_files = TRUE){
  if(write_files == TRUE){
    annotate_variants_with_intermediates(bed_dir_path = bed_dir_path, variant_path = variant_path, cores = cores)
  }else{
    annotate_variants(bed_dir_path = bed_dir_path, variant_path = variant_path, cores = cores)
  }
}

## Functions for annotating from vector of file paths ----------------
##' Write intermediate cols and horizontally concatenate output.
##'
##' @family annotation functions WITH intermediate file output
##' @param bed_dir_path A character.
##' @param variant_path A character.
##' @return NULL.
##' @export
