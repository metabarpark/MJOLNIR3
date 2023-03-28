# LOKI: LULU Overseeing with Kinship Identification

# This function is a convenient wrapper of LULU for the MJOLNIR metabarcoding pipeline.
# It starts from the combined dataset of abundances and taxonomy from the previous step (FRIGGA): LIBR.All_MOTUs.tsv
# Then a match list of representative MOTU sequences is created using VSEARCH and saved as a txt file.
# Then MOTUs which are potential pseudogenes are labelled and removed using LULU.
# The output includes the LIBR.match_list.txt file and 3 CSV files:
# - A final curated metabarcoding dataset: LIBR.Curated_LULU.tsv
# - A dataset of discarded MOTUs: LIBR.Discarded_LULU.tsv
# - A file with informaton on the fate of discarded MOTUs (with IDs of putative mother sequences): LIBR.Deleted_LULU_fate.tsv

mjolnir7_LOKI <- function(lib,min_id = .84, vsearchpath = NULL){

  # lib is the name of the library to be processed. Usually a four-character uppercase name.
  # Input file name must be in the format: LIBR.All_MOTUs.tsv. Then lib must be = "LIBR"
  # min_id is the minimum identity between two sequences to be kept in the match_list output. Default: 0.84

  if (is.null(vsearchpath)) vsearchpath <- "~/vsearch-2.22.1/bin/"
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(old_path, path.expand(vsearchpath), sep = ":"))

  message("LOKI will produce a pairwise match list for LULU.")
  system(paste0("cat ",lib,"_ODIN_part_??.fasta > ",lib,"_LOKI.fasta"))
  system(paste0("vsearch --usearch_global ",lib,"_LOKI.fasta --db ",lib,"_LOKI.fasta --self --id ",min_id," --iddef 1 --userout ",lib,"_LOKI_match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"),intern=T,wait=T)
  message("LOKI will now remove the pseudogenes with LULU.")

  suppressPackageStartupMessages(library("lulu"))
  suppressPackageStartupMessages(library("dplyr"))

  #Load the dataset
  db <- read.table(paste0(lib,"_FRIGGA.tsv"),sep="\t",head=T,stringsAsFactors = F)
    # Select sample abundance columns
  sample_cols <- grep("sample",names(db))
  otutable_name <-db[,sample_cols]
  rownames(otutable_name) <- db$id

  #Load the matchlist
  matchlist_name <- read.csv(paste0(lib,"_LOKI_match_list.txt"),sep="\t",head=F,stringsAsFactors = F)

  #Run LULU
  curated_result <- lulu(otutable_name, matchlist_name)

  #Get discarded table:
  discarded_db <- db[db$id %in% curated_result$discarded_otus,]
  otus_discarded <- curated_result$discarded_otus
  num_discarded <- nrow(discarded_db)
  write.table(discarded_db,paste0(lib,"_LOKI_Discarded.tsv"),row.names = F,sep="\t",quote = F)

  #Get curated table:
  curated_db <- db[db$id %in% curated_result$curated_otus,]
  curated_db <- curated_db[order(curated_db$id),]
  curated_db[,sample_cols] <- curated_result$curated_table
  curated_db$total_reads <- rowSums(curated_result$curated_table)
  write.table(curated_db,paste0(lib,"_LOKI_Curated.tsv"),row.names = F,sep="\t",quote = F)

  #Get fate of deleted taxa
  deleted_otu_fate <- (curated_result$otu_map[curated_result$otu_map$curated=="merged",])
  deleted_otu_fate$original <- ""
  deleted_otu_fate$id_removed <- rownames(deleted_otu_fate)
  for (i in 1:nrow(deleted_otu_fate)){
    deleted_otu_fate$original[i] <- db$SCIENTIFIC_NAME[db$id==rownames(deleted_otu_fate)[i]]
  }
  parents_info <- db[db$id%in%deleted_otu_fate$parent_id,c("id","SCIENTIFIC_NAME","superkingdom_name","kingdom_name","phylum_name","class_name","order_name","family_name")]
  names(parents_info)[1:2] <- c("parent_id","parent_taxo")
  deleted_otu_fate <- left_join(deleted_otu_fate,parents_info,by="parent_id")
  write.table(deleted_otu_fate,paste0(lib,"_LOKI_Deleted_fate.tsv"),row.names = F,sep="\t",quote = F)

  save(file = "summary_LOKI.RData",list = c("num_discarded","otus_discarded"))

  message("LOKI is done. He kept ",nrow(curated_db)," MOTUs in the curated database, which He stored in file: ",paste0(lib,"_LOKI_Curated.tsv"))
  message("LOKI discarded ",nrow(discarded_db)," MOTUs, which He stored in file: ",paste0(lib,"_LOKI_Discarded.tsv"))
  message("LOKI stored the fate of the discarded MOTUs in file: ",paste0(lib,"_LOKI_Discarded.tsv"))
}
