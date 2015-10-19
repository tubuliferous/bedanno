library(data.table)
library(colorout)
library(stringr)


move_file <- function(path){
  to_path <- "/Users/tubuliferous/Desktop/ENCODE_dump/K562/"
  system(paste("mv", path, paste(to_path, basename(path), sep="")))  
}


manifest <- fread("/Users/tubuliferous/Dropbox/eQTL_publication_work/Annotations/ENCODE_dump/metadata.tsv")
anno_dir <- "/Users/tubuliferous/Desktop/ENCODE_dump/ENCODE_bed_files/"
file_paths <- list.files(anno_dir)

names(manifest) <- str_replace_all(names(manifest), " ", "_")
subset(manifest, Biosample_term_name=="K562")[, "File_accession", with=F]

k562_accesssions <- subset(manifest, Biosample_term_name=="K562")$File_accession
k562_accesssions_narrow <- paste(k562_accesssions, ".narrowPeak.gz", sep = "")
k562_accesssions_broad <- paste(k562_accesssions, ".broadPeak.gz", sep = "")

full_file_paths <- paste(anno_dir, file_paths, sep="")
k562_file_paths_narrow <- full_file_paths[file_paths %in% k562_accesssions_narrow]
k562_file_paths_broad <- full_file_paths[file_paths %in% k562_accesssions_broad]


sapply(k562_file_paths_broad, move_file)


subset(manifest, Biosample_term_name=="K562", )


str_dete$t(manifest$Experiment_target, "POL2")
manifest[manifest$Biosample_term_name=="K562" & str_detect("RAD1", manifest$Experiment_target), ]

sub_k562 <- subset(manifest, Biosample_term_name=="K562")

sub_k562[str_detect(sub_k562$Experiment_target, "RAD"), ]


hepg2 <- manifest[manifest$"Biosample term name" == "HEPG2", ]