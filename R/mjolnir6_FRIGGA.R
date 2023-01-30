# FRIGGA: Final Recount and Integration of Generated Genealogies and Abundances

## Script for combining an abundance TSV file from ODIN with a taxonomy-annotated TSV file from THOR.
## The script will have the name of the library (typically 4 characters, e.g. LIBR) as input.
## The taxonomy file must be called LIBR.ecotag.fasta.annotated.tsv
## The abundance file must be called LIBR.SWARM_output.counts.tsv
## The separator characters for both the taxonomy and abundance tsv files must be tabs.
## The abundances file must include sample columns names starting with "sample."
## The output file is a TSV file called LIBR.All_MOTUs.tsv
## FRIGGA deprecates the owi_combine function from previous pipelines (Project Metabarpark, 2016).
## By Owen S. Wangensteen

mjolnir6_FRIGGA <- function(lib=NULL){

  message("FRYGGA will produce a combined file.")

  infile=paste0(lib,"_THOR_annotated.tsv")
  abundances=paste0(lib,"_ODIN_counts.tsv")
  outfile=paste0(lib,"_FRIGGA.tsv")
  if (!file.exists(abundances)) {
    abundances=paste0(lib,"_ODIN_ESV.tsv")
  }

  message("FRYGGA is reading the ecotag-annotated database from THOR...")
  ecotag_db <- read.table(infile,sep="\t",head=T,stringsAsFactors=F)
  message("FRYGGA has read the taxonomy database, with ", nrow(ecotag_db)," total MOTUs.")
  # Delete "None" from the taxonomy database
  ecotag_db[ecotag_db=="None"] <- ""

  message("FRYGGA is reading the abundances database from ODIN...")
  abun_db <- read.table(abundances,sep="\t",head=T,stringsAsFactors=F)
  n_samples <- length(grep("sample",names(abun_db)))
  message("FRYGGA has read the abundances database, including ", nrow(abun_db)," total MOTUs and ",n_samples," samples.")

  # Merge databases
  names(abun_db)[names(abun_db)=="ID"] <- "id"
  names(abun_db)[names(abun_db)=="NUC_SEQ"] <- "sequence"

  db <- merge(ecotag_db,abun_db,by="id")

  names(db)[grep("sample",names(db))] <- gsub("MERGED_sample.","",names(db)[grep("sample",names(db))])
  db$COUNT <- rowSums(db[,grep("sample",names(db))])

  write.table(db,outfile,sep="\t",quote=F,row.names=F)
  message("FRYGGA is done. File ", outfile, " written, including ",nrow(db)," MOTUs with ",sum(db$total_reads)," total reads in ",n_samples," samples.")
  message("(",sum(db$COUNT>1)," non-singletons MOTUs).")
}

