#' @docType package
#' @name bedanno
#' @title Annotate BED-like variant files.
#' @import dplyr
NULL
#' @import data.table
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to bedanno! Both variants and annotations files should follow the BED schema for the first three columns.")
}

# Utility Functions ---------------------------------------------

#' Read gzipped file into data.table.
#'
#' @family utility functions
#' @param path A character.
#' @param sep A character.
#' @return data.table
gzfread <- function(path, sep, out_dir = NULL){
  if(!stringr::str_detect(path, ".gz$")){
    data.table::fread(path, sep = sep)
  }else{
    if(is.null(out_dir)){
      gunzipped_path <- stringr::str_replace(path, ".gz$", "")
    }else{
      file_name <- stringr::str_replace(basename(path, ".gz$", ""))
      gunzipped_path <- paste(out_dir, "/", filename, sep="")
    }

    gunzipped_path <- stringr::str_replace(path, ".gz$", "")
    if(file.exists(gunzipped_path)){
      data.table::fread(gunzipped_path, sep = "")
    }else{
      gunzipped <- data.table::fread(R.utils::gunzip(path, remove=F), sep = sep)
    }
  }
}

#' Import BED file.
#'
#' @family utility functions
#' @param bed_path A character.
#' @return data.table
get_bed <- function(bed_path){
  paste("Reading ", basename(bed_path), sep="") %>% print
  this_bed <- read_delim(bed_path, delim="\t", col_names=FALSE)
  # Alternative: this_bed <- gzfread(bed_path, sep="\t")
  this_bed <- data.table::data.table(this_bed)
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
  file_list <- paste(dir_path, list.files(dir_path), sep="")
  return(file_list)
}

# Main Functions -----------------------------------------------

#' Generate binary annotation column.
#'
#' @family annotation functions with NO intermediate file output
#' @param bed_file_path A character.
#' @param variant_table A data.table.
#' @return Annotation column.
#' @export
get_anno_col <- function(bed_file_path, variant_table){
  bed      <- get_bed(bed_file_path)
  data.table::setkey(bed, CHROM, START, STOP)
  overlap_table <- foverlaps(variant_table, bed, type="within")
  distinct_overlap_table <- unique(overlap_table, by=c("CHROM", "i.START", "i.STOP"))
  anno_column <- rep.int(0, nrow(distinct_overlap_table)) # Use rep.int instead of rep for speed
  anno_column[!is.na(distinct_overlap_table$START)] <- 1
  anno_column %>% return
}

#' Annotate a data.table with binary annotations.
#'
#' @family annotation functions with NO intermediate file output
#' @param bed_dir_path A character.
#' @param variant_path A data.table.
#' @return Annotated data.frame.
#' @export
annotate_variants <- function(bed_dir_path, variant_path, cores=1){
  elapsed_get_variant_table <- system.time(variant_table <- get_variant_table(variant_path))
  print(paste("Variant table read in", elapsed_get_variant_table["elapsed"], "seconds."))
  bed_paths <- get_file_paths(bed_dir_path)
  elapsed_anno_list <- (system.time(anno_list <- parallel::mclapply(bed_paths, get_anno_col, variant_table, mc.cores=cores)))
  print(paste("Annotation matrix generated in", elapsed_anno_list["elapsed"], "seconds"))
  names(anno_list) <- basename(bed_paths)
  return(data.frame(variant_table, anno_list))
}

# -------------------------------------------

#' Write intermediate annotation column.
#'
#' @family annotation functions WITH intermediate file output
#' @param bed_file_path A character.
#' @param variant_table A data.table.
#' @param anno_path A character.
#' @return NULL.
#' @export
write_anno_col <- function(bed_file_path, variant_table, anno_path){
  bed      <- get_bed(bed_file_path)
  data.table::setkey(bed, CHROM, START, STOP)
  overlap_table <- foverlaps(variant_table, bed, type="within")
  distinct_overlap_table <- unique(overlap_table, by=c("CHROM", "i.START", "i.STOP"))
  anno_column <- rep.int(0, nrow(distinct_overlap_table)) # Use rep.int instead of rep for speed
  anno_column[!is.na(distinct_overlap_table$START)] <- 1
  anno_column <- data.frame(anno_column)
  names(anno_column) <- basename(bed_file_path)
  write.table(anno_column, file = paste(anno_path, "/", basename(bed_file_path), ".anno", sep=""), col.names = TRUE, quote = FALSE, row.names = FALSE)
}

#' Write intermediate annotation columns from dir.
#'
#' @family annotation functions WITH intermediate file output
#' @param bed_dir_path A character.
#' @param variant_path A character.
#' @return NULL.
write_variant_annotations <- function(bed_dir_path, variant_path, cores = 1){
  anno_dir_path <- file.path(dirname(variant_path), paste(basename(variant_path), "annos", sep="."))
  dir.create(anno_dir_path, showWarnings = FALSE)
  elapsed_get_variant_table <- system.time(variant_table <- get_variant_table(variant_path))
  write.table(variant_table, paste(variant_path, ".reordered", sep=""), col.names = TRUE, quote = FALSE, row.names = FALSE)
  print(paste("Variant table read in", elapsed_get_variant_table["elapsed"], "seconds."))
  bed_paths <- get_file_paths(bed_dir_path)
  elapsed_anno <- (system.time(parallel::mclapply(bed_paths, write_anno_col, variant_table, anno_dir_path, mc.cores = cores)))
  print(paste("Annotation files written in", elapsed_anno["elapsed"], "seconds"))
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
  # Use UNIX paste command
  sys_command <- paste("paste ", paths, " > ", variant_path, ".annotated", sep="")
  print("Horizontal concatentation system command:")
  print(sys_command)
  system(sys_command)
}

#' Write intermediate cols and horizontally concatenate output.
#'
#' @family annotation functions WITH intermediate file output
#' @param bed_dir_path A character.
#' @param variant_path A character.
#' @return NULL.
#' @export
annotate_variants_with_intermediates <- function(bed_dir_path, variant_path, cores = 1){
  write_variant_annotations(bed_dir_path = bed_dir_path, variant_path = variant_path, cores = cores)
  anno_col_dir_path <- file.path(dirname(variant_path), paste(basename(variant_path), "annos", sep="."))
  formatted_var_path <- paste(variant_path, ".reordered", sep="") #From

  dir.create(paste(dirname(bed_dir_path), "/", "unzipped_bed_dir", sep = ""), showWarnings = FALSE)
  bed_dir_path 

  horizontal_concat_annos(variant_path = formatted_var_path, anno_col_dir_path = anno_col_dir_path)
}



