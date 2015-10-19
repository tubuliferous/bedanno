require(colorout)
library(dplyr)
library(devtools)
install_github("tubuliferous/bedanno")
library(bedanno)
ls("package:devtools")
ls("package:bedanno")
library(data.table)

meta <- fread("/Users/tubuliferous/Dropbox/mygithub/bedanno/inst:extdata/encode_dump/metadata.tsv")
encode_files <- list.files("/Users/tubuliferous/Dropbox/mygithub/bedanno/inst:extdata/encode_dump/anno_encode_files/")

anno_vars(bed_dir_path = "/Users/tubuliferous/Desktop/bedanno_test/elements_folder/", variant_path = "/Users/tubuliferous/Desktop/bedanno_test/variants", write_files = FALSE)
