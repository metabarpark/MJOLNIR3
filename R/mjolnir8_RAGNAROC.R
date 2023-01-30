# RAGNAROC: Replace AGnomens with Names And Recover Original Codification
# Sort_MOTUs options available: "id" (default), "taxonomy", "abundance"
# To remove bacterial MOTUs use remove_bacteria = T
# To remove contaminant MOTUs use remove_contamination = T.
# When remove_contamination = T a list of contaminants must be provided in contamination_file (a text file with one line for each contaminant scientific_name)


mjolnir8_RAGNAROC <- function(lib,metadata_table="",output_file="",output_file_ESV="",min_reads=2,min_relative=1/50000,
                              remove_bacteria=F,remove_contamination=F,contamination_file="contaminants.txt",
                              ESV_within_MOTU=F,blank_col="BLANK",blank_tag=T,remove_numts=F,cores=1,blank_relative=0.1){

  message("RAGNAROC is coming. Original sample names will be recovered.")
  suppressPackageStartupMessages(library(tidyr))
  if (output_file == "") {
    output_file <- paste0(lib,"_RAGNAROC_final_dataset.tsv")
  }
  if (output_file_ESV == "") {
    output_file_ESV <- paste0(lib,"_RAGNAROC_final_dataset_ESV.tsv")
  }

  if (file.exists("summary_FREYJA.RData")) {
    load("summary_FREYJA.RData")
    rownames(variables_FREYJA) <- variables_FREYJA$variable
    FREYJA <- T
  } else {
    FREYJA <- F
  }
  if (file.exists("summary_HELA.RData")) {
    load("summary_HELA.RData")
    HELA <- T
  } else {
    HELA <- F
  }
  if (file.exists("summary_ODIN.RData")) {
    load("summary_ODIN.RData")
    ODIN <- T
  } else {
    ODIN <- F
  }
  if (file.exists("summary_LOKI.RData")) {
    load("summary_LOKI.RData")
    LOKI <- T
  } else {
    LOKI <- F
  }

  # Load the data set
  # two data sets can be used:
  # MOTU/ESV + taxa info with LULU/without LULU
  # ESV within MOTU only available for ODIN_algorithm DnoisE_SWARM/SWARM_DnoisE

  if (LOKI | file.exists(paste0(lib,"_LOKI_Curated.tsv"))) {
    db <- read.csv(paste0(lib,"_LOKI_Curated.tsv"), sep="\t",head=T,stringsAsFactors = F)
  } else {
    db <- read.csv(paste0(lib,"_FRIGGA.tsv"), sep="\t",head=T,stringsAsFactors = F)
  }
  if (ESV_within_MOTU) {
    ESV_data_initial <- read.csv(paste0(lib,"_ODIN_ESV.tsv"), sep="\t",head=T,stringsAsFactors = F)
  }

  # Remove bacteria
  if (remove_bacteria) {
    message("RAGNAROC is removing bacterial MOTUs now.")
    bacteria_removed <- sum(db$superkingdom_name == "Prokaryota" | db$SCIENTIFIC_NAME == "root")
    db <- db[(db$superkingdom_name != "Prokaryota" & db$SCIENTIFIC_NAME != "root"),]
  }

  # Remove contamination
  if (remove_contamination){
    message("RAGNAROC is removing contaminant MOTUs now.")
    contamination <- readLines(contamination_file)
    db <- db[!((db$SCIENTIFIC_NAME %in% contamination) |
                         (db$phylum_name %in% contamination) |
                         (db$class_name %in% contamination) |
                         (db$order_name %in% contamination) |
                         (db$family_name %in% contamination) |
                         (db$genus_name %in% contamination)) ,]
  }

  # Load the metadata_table
  if (metadata_table=="") metadata_table <- paste0(lib,"_metadata.tsv")
  sample_db <- read.table(metadata_table,sep="\t",head=T,stringsAsFactors = F)

  message("data loaded")

  # if ESV_within_MOTU remove the MOTUs that LOKI has removed and bacteria filter
  if (ESV_within_MOTU) {
    ESV_data_initial <- ESV_data_initial[ESV_data_initial$MOTU%in%db$id,]
  }

  # Select sample abundance columns
  sample_cols <- grep("sample",names(db))
  initial_no_sample_cols <- length(sample_cols)
  sample_names <- names(db[sample_cols])

  # Change agnomens by original names
  new_sample_names <- sample_db$original_samples[match(sample_names,sample_db$mjolnir_agnomens)]
  new_sample_names[is.na(new_sample_names)] <- gsub("^","EMPTY",as.character(c(1:sum(is.na(new_sample_names)))))
  colnames(db)[sample_cols] <- new_sample_names

  # get negatives/blanks
  if (!ESV_within_MOTU) {
    neg_samples <- db[,sample_cols[grepl(paste0(sample_db$original_samples[as.character(sample_db[,blank_col])==as.character(blank_tag)],collapse = "|"),new_sample_names)]]
  }

  # remove negs and empties
  db <- db[,!grepl(paste(c(negatives,blank,"EMPTY"),collapse = "|"),colnames(db))]

  # correct sample identifiers
  sample_names <- new_sample_names[!grepl(paste(c(negatives,blank,"EMPTY"),collapse = "|"),new_sample_names)]
  sample_cols <- match(sample_names,colnames(db))

  # same for ESVs
  if (ESV_within_MOTU) {
    # Select sample abundance columns
    sample_cols_ESV <- grep("sample",names(ESV_data_initial))
    sample_names_ESV <- gsub("MERGED_sample.","",names(ESV_data_initial[sample_cols_ESV]))

    # Change agnomens by original names
    new_sample_names_ESV <- sample_db$original_samples[match(sample_names_ESV,sample_db$mjolnir_agnomens)]
    new_sample_names_ESV[is.na(new_sample_names_ESV)] <- gsub("^","EMPTY",as.character(c(1:sum(is.na(new_sample_names_ESV)))))
    names(ESV_data_initial)[sample_cols_ESV] <- new_sample_names_ESV

    neg_samples <- ESV_data_initial[,sample_cols_ESV[grepl(paste0(sample_db$original_samples[as.character(sample_db[,blank_col])==as.character(blank_tag)],collapse = "|"),new_sample_names_ESV)]]

    # remove neg samples
    if (dim(neg_samples)[2]>0) {
      # remove negs and empties
      ESV_data_initial <- ESV_data_initial[,!grepl(paste(c(negatives,blank,"EMPTY"),collapse = "|"),colnames(ESV_data_initial))]
      # correct sample identifiers
      sample_names_ESV <- new_sample_names[!grepl(paste(c(negatives,blank,"EMPTY"),collapse = "|"),new_sample_names)]
      sample_cols_ESV <- match(sample_names_ESV,colnames(ESV_data_initial))
    }
  }

  # db[1,"09_12_M2_A_C"] <- 10000000000000



  # Filter 1. remove any MOTU for which abundance in the blank or negative controls was higher than 10% of its total read abundance
  # remove blank and NEG samples
  if (dim(neg_samples)[2]>0) {
    message("RAGNAROC will remove any MOTU for which abundance in the blank or negative controls was higher than 10% of its total read abundance")
    neg_reads <- rowSums(neg_samples)
    if (!ESV_within_MOTU) {
      sample_reads <- rowSums(db[,sample_cols])
      data_neg_filt_deleted <- db[neg_reads/(sample_reads+neg_reads) > blank_relative,]
      db <- db[!neg_reads/(sample_reads+neg_reads) > blank_relative,]
    } else {
      sample_reads <- rowSums(ESV_data_initial[,sample_cols_ESV])
      data_neg_filt_deleted <- ESV_data_initial[neg_reads/(sample_reads+neg_reads) > blank_relative,]
      ESV_data_initial <- ESV_data_initial[!neg_reads/(sample_reads+neg_reads) > blank_relative,]
      db <- db[db$id %in% unique(ESV_data_initial$MOTU),]
    }
    message("Blank correction finished")
  }

  # Filter 2. Apply a minimum relative abundance threshold for each sample, setting to zero any abundance below min_relative of the total reads of this sample
  # it also applies a min_reads filter
  message("RAGNAROC is applying a relative abundance filter. MOTUs with less than ",min_relative," relative abundance will be removed from each sample.")

  relabund <- function(x,min_relative) if (sum(x)>0) x/sum(x) < min_relative else FALSE
  if (!ESV_within_MOTU) {
    rownames(db) <- db$id

    change_matrix <- do.call("cbind",apply(db[,sample_cols], 2, relabund, min_relative=min_relative)) & db[,sample_cols]>0

    relabund_changed <- data.frame(id_modified = rownames(change_matrix[rowSums(change_matrix)>0,]),
                                   samples = vapply(rownames(change_matrix[rowSums(change_matrix)>0,]), function(x,change_matrix){
                                     return(paste(colnames(change_matrix)[change_matrix[rownames(change_matrix)==x,]]))
                                   }, FUN.VALUE = "string", change_matrix = change_matrix))
    db[,sample_cols][change_matrix] <- 0
    db$COUNT <- rowSums(db[,sample_cols])
    message("Filter of minimum relative abundance finished")
    message("RAGNAROC is removing MOTUs with less than ",min_reads," total reads.")
    db <- db[db$COUNT >= min_reads,]
  } else {
    rownames(ESV_data_initial) <- ESV_data_initial$ID

    change_matrix <- do.call("cbind",apply(ESV_data_initial[,sample_cols_ESV], 2, relabund, min_relative=min_relative)) & ESV_data_initial[,sample_cols_ESV]>0

    relabund_changed <- data.frame(ESV_id_modified = rownames(change_matrix[rowSums(change_matrix)>0,]),
                                   samples = vapply(rownames(change_matrix[rowSums(change_matrix)>0,]), function(x,change_matrix){
                                     return(paste(colnames(change_matrix)[change_matrix[rownames(change_matrix)==x,]],collapse = "|"))
                                   }, FUN.VALUE = "string", change_matrix = change_matrix))
    ESV_data_initial[,sample_cols_ESV][change_matrix] <- 0

    ESV_data_initial$COUNT <- rowSums(ESV_data_initial[,sample_cols_ESV])

    message("RAGNAROC is removing MOTUs with less than ",min_reads," total reads.")
    ESV_data_initial <- ESV_data_initial[ESV_data_initial$COUNT >= min_reads,]

    # remove numts
    if (remove_numts) {
      message("numts will be removed")
      no_ESV_before_numts <- dim(ESV_data_initial)[1]
      lengths <- nchar(as.vector(ESV_data_initial$NUC_SEQ))
      ESV_data_initial <- ESV_data_initial[(lengths-313)%%3 == 0,]
      lengths <- nchar(as.vector(ESV_data_initial$NUC_SEQ))

      no_numts_data <- c()
      numts_seqs <- c()

      number_of_motus <- length(unique(ESV_data_initial$MOTU))
      motu_taxa <- data.frame("id" = db$id, "Metazoa" = c(db$kingdom_name == "Metazoa" & !is.na(db$kingdom_name)))
      numts_ESV <- parallel::mclapply(1:number_of_motus,function(i,ESV_data_initial,motu_taxa){
        motu <- unique(ESV_data_initial$MOTU)[i]
        datas <- ESV_data_initial[ESV_data_initial$MOTU==motu,]
        is_metazoa <- motu_taxa$Metazoa[motu_taxa$id==as.character(motu)]
        datas_length <- nchar(as.vector(datas$NUC_SEQ))
        newlist <- numts(datas, is_metazoa = is_metazoa, motu = motu, datas_length = datas_length)
        return(newlist)
      },ESV_data_initial=ESV_data_initial,motu_taxa=motu_taxa,mc.cores = cores)
      numts_ESV <- do.call("rbind",numts_ESV)
      ESV_data_initial <- ESV_data_initial[ESV_data_initial$ID %in% numts_ESV$id,]
      message("numts removed")
    }
    db <- db[db$id %in% unique(ESV_data_initial$MOTU),]
    # compute new abundances of MOTUs from ESV
    motu_abund <- parallel::mclapply(db$id, FUN = function(x,ESV_data_initial,sample_cols_ESV){
      data_motu <- ESV_data_initial[as.character(ESV_data_initial$MOTU) == x,sample_cols_ESV]
      if (dim(data_motu)[1]>1) {
        return(colSums(data_motu))
      } else {
        return(data_motu)
      }
    },ESV_data_initial=ESV_data_initial,sample_cols_ESV=sample_cols_ESV,mc.cores = cores)
    db[,sample_cols] <- do.call("rbind",motu_abund)
  }

  # Write final table
  write.table(db,output_file,row.names = F,sep="\t",quote = F)
  if (ESV_within_MOTU){
    write.table(ESV_data_initial,output_file_ESV,row.names = F,sep="\t",quote = F)
  }
  message("After RAGNAROC, MJOLNIR is done. File: ",output_file, " written with ",nrow(db_new), " MOTUs and ",sum(db_new$total_reads)," total reads.")

  sink("RAGNAROC_summary_report.txt")
  cat(paste0("Dear friend,\n",
             "you have succesfully arrived at the end of RAGNAROC. You've meet gods and took their help to twist the data to your will.\n",
             "After RAGNAROC the rest is up to you. Don't lose the faith in your experiment, the end is near but new paths will open below your feet.\n",
             "Please don't forget to cite and thank the two dwarfs Cindri and Brok, AKA Owen and Adria, for the forge of my self.\n",
             "MJOLNIR.\n",
             "P.S.: See below for a small summary of your journey.\n"))
  if (FREYJA) {
    if (as.logical(variables_FREYJA["demultiplexed",2])) {
      cat(paste0("You started FREYJA with your samples allready demultiplexed and with the following sequences for each file \n"))
      do.call("rbind",before_FREYJA)
    } else{
      cat(paste0("You started FREYJA with the following sequences for each file \n"))
      do.call("rbind",lapply(before_FREYJA,function(x)do.call("rbind",x)))
    }
    cat(paste0("You used ",variables_FREYJA["cores",2]," cores to aling your sequences. You choosed those sequences with a quality score of more than ",variables_FREYJA["score_obialign",2],".\n",
               "You assign each sequence to a sample name and removed the primer's sequences.\n",
               "Finally in FREYJA you just kept those sequences with A, G, T or C's and with a sequence length between ",variables_FREYJA["Lmin",2]," and ",variables_FREYJA["Lmax",2]," bp.\n",
               "The resulting files had the following stats:\n"))
    as.data.frame(pivot_wider(do.call("rbind",after_FREYJA),names_from = "version",values_from = "num_seqs"))
  } else {
    cat("Sorry but I couldn't find a summary of your FREYJA process")
  }
  if (HELA) {
    cat(paste0("HELA removed the chimeras with the uchime_denovo algorithm and kept for each sample the following number of non-chimeras:\n"))
    after_HELA
  } else {
    cat("Sorry but I couldn't find a summary of your HELA process")
  }
  if (ODIN) {
    cat(paste0("ODIN was used to obtain meaningful units. In your case you chose the ",algorithm," algorithm.\n"))
    if (algorithm=="dnoise_swarm" | algorithm=="dnoise") {
      cat(paste0("ODIN used DnoisE to obtain the ESV's of your samples running within them with the following options:\n"))
      if (COI) {
        cat(paste0("Entropy correction with sequences delimited to a multiple of 313bp, alpha ",alpha," and minimum number of reads of ",min_reads_ESV,"\n"))
      } else {
        cat(paste0("Alpha ",alpha," and minimum number of reads of ",min_reads_ESV,"\n"))
      }
    }
    if (algorithm=="dnoise_swarm"  | algorithm=="swarm" | algorithm=="swarm_dnoise") {
      cat(paste0("ODIN joined all the sequences, obtained the unique ones and applied swarm to obtain the MOTUs. Before SWARM you had",after_2_ODIN$values[after_ODIN$version==seq_id]," sequences and at the end you obtained the following stats: \n"))
      after_4a_ODIN
    } else{
      cat(paste0("The samples were then grouped and the unique sequences obtained being ",after_2_ODIN$values[after_ODIN$version==seq_id]," sequences in total.\n"))
    }
    if (algorithm=="swarm_dnoise") {
      cat(paste0("ODIN used DnoisE to obtain the ESV's of your samples running within them with the following options:\n"))
      if (COI) {
        cat(paste0("Entropy correction with sequences delimited to a multiple of 313bp, alpha ",alpha," and minimum number of reads of ",min_reads_ESV,"\n"))
      } else {
        cat(paste0("Alpha ",alpha," and minimum number of reads of ",min_reads_ESV,"\n"))
      }
    }
  } else {
    cat("Sorry but I couldn't find a summary of your ODIN process")
  }
  if (LOKI) {
    cat(paste0("LOKI used LULU to search for potential pseudogenes and found ",num_discarded," OTUs that were discarded.\n"))
    after_HELA
  } else {
    cat("Sorry but I couldn't find a summary of your LOKI process")
  }
  cat(paste0("During RAGNAROC some filters were applied."))
  if (remove_bacteria) cat(paste0(bacteria_removed," bacteria were removed"))
  if (remove_contamination) cat(paste0("contaminations were removed"))
  if (dim(neg_samples)[2]>0) cat(paste0(data_neg_filt_deleted,ifelse(ESV_within_MOTU," ESV"," MOTU")," were removed by neg/blank filter"))
  cat(paste0("The relative abundance filter of ",min_relative," within samples had effect on the following id's and samples:"))
  relabund_changed
  if (ESV_within_MOTU&remove_numts) cat(paste0("The numts filter found ",numts_ESV," numts that were removed."))
  sink()
}

compare.DNA <- function(x,y){
  as.integer(x) == as.integer(y)
}
numts<-function(datas, is_metazoa=FALSE, motu, datas_length)
{
  library(Biostrings)
  library(stringr)

  # compare only mitochondrial genetic code
  mitochondrial_GC <- c(2,3,4,5,7,11,12,14,15,16,17,18)
  # START
  motu_name = motu
  datas$NUC_SEQ<-as.character(datas$NUC_SEQ)

  # remove sequences with different length than the seed
  if (sum(datas$ID==motu)==0) { # if the seed has been deleted in previous steps take the first more abundant
    motu = datas$ID[which(datas$COUNT==max(datas$COUNT,na.rm = TRUE))[1]]
  }
  correct_length <- datas_length[datas$ID==motu]
  datas <- datas[datas_length==correct_length,]

  # remove misaligned sequences (more than 30 differences between a sequence and
  # the seed)
  misaligned_seqs <- c()
  motu_seq <- DNAString(datas$NUC_SEQ[datas$ID == motu])
  for (i in 1:dim(datas)[1]) {
    if(sum(!compare.DNA(motu_seq,DNAString(datas$NUC_SEQ[i])))>=30){
      misaligned_seqs <- c(misaligned_seqs,i)
    }
  }
  if (length(misaligned_seqs)>0) {
    datas <- datas[-misaligned_seqs,]
  }

  # look for the best genetic code, this is the one with less stop codons in all sequences
  # the number of codons stop is multiplied by the number of count of the sequence.
  stops<-matrix(NA,dim(datas)[1],20)
  aa_xung<-matrix(NA,dim(datas)[1],20)
  seq<-DNAStringSet(datas$NUC_SEQ)
  seq<-DNAStringSet(seq,start=2,end=nchar(seq[1]))

  for (qq in mitochondrial_GC){
    code<-getGeneticCode(as.character(GENETIC_CODE_TABLE$id[qq]))
    trans<-translate(seq,genetic.code=code)

    # for (k in 1:dim(datas)[1]){
    # nstops <- stringr::str_count(as.character(trans[k]),fixed("*"))
    # stops[k,qq] <- nstops * datas$COUNT[k]
    # }
    nstops <- apply(data.frame(as.character(trans)), 1, function(x){stringr::str_count(x,fixed("*"))})
    stops[,qq] <- nstops * datas$COUNT

  }

  goodcodes<-which(colSums(stops)==min(colSums(stops),na.rm = T))

  # if more than one code have been chosen as good code choose the first as the best.
  # However, if the MOTU is a Metazoan and has 313bp length we check the 5 well preserved aa
  # and the code with less changes per read is the chosen. Also remove sequences
  # with changes in thees positions

  if (is_metazoa & (correct_length == 313)){

    for (qq in 1:length(goodcodes))
    {
      code<-getGeneticCode(as.character(GENETIC_CODE_TABLE$id[goodcodes[qq]]))
      trans<-translate(seq,genetic.code=code)
      for (k in 1:dim(datas)[1])
      {
        aa<-strsplit(as.character(trans[k]),split="")
        aa<-unlist(aa)
        bad_aa<-0
        if (aa[20]!="H") {bad_aa<-bad_aa+datas$COUNT[k]} # the number of errors counted are the same as the number of counts of the seq.
        if (aa[23]!="G") {bad_aa<-bad_aa+datas$COUNT[k]}
        if (aa[32]!="N") {bad_aa<-bad_aa+datas$COUNT[k]}
        if (aa[81]!="D") {bad_aa<-bad_aa+datas$COUNT[k]}
        if (aa[95]!="G") {bad_aa<-bad_aa+datas$COUNT[k]}
        aa_xung[k,goodcodes[qq]]<-bad_aa
      }
    }

    goodcodes <- goodcodes[which(colSums(aa_xung)[goodcodes]==min(colSums(aa_xung)[goodcodes],na.rm = T))]
    bestcode <- goodcodes[1]
    bestcodename <- GENETIC_CODE_TABLE$name[bestcode]
    goodcodesnames <- GENETIC_CODE_TABLE$name[goodcodes]
    flag <- stops[,bestcode]>0 | aa_xung[,bestcode]>0
  } else {
    bestcode<-which(colSums(stops)==min(colSums(stops),na.rm = T))[1]
    bestcodename<-GENETIC_CODE_TABLE$name[bestcode]
    goodcodesnames <- GENETIC_CODE_TABLE$name[goodcodes]
    flag <- stops[,bestcode]>0
  }

  # numts
  if (sum(flag)>0) {
    numts_seqs <- data.frame("motu" = motu_name, "id" = datas$ID[flag],
                             "genetic_code" = bestcodename,
                             "similar_codes" = paste(goodcodesnames, collapse = " | "))
  } else {
    numts_seqs <- c()
  }

  # # remove numts
  # datas <- datas[(flag==FALSE),]
  #
  # newlist <- list("no_numts_data" = datas, "numts_seqs" = numts_seqs)

  return(numts_seqs)
}
