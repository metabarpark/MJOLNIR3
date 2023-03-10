# HELA: Hierarchical Elimination of Lurking Artifacts

# This function uses the uchime_denovo algorithm implemented in VSEARCH to remove chimaeric sequences from the dataset.
# HELA works in a sample-by-sample basis. HELA will process all individual fasta files in the current folder matching the pattern XXXX_sample_XXX.fasta.
# This allows for parallel computing, significantly decreasing calculation times.
# HELA can optionally remove singleton sequences (default: remove_singletons=T), so that the computing time for clustering with ODIN will be considerably reduced.
# This is possibly a good idea for very large datasets (with > 5 million unique sequences before clustering)
# The final dataset output is in VSEARCH format, so it can be directly fed into SWARM (ODIN).

mjolnir3_HELA <- function(lib, cores, vsearchpath = NULL){

  suppressPackageStartupMessages(library(parallel))
  old_path <- Sys.getenv("PATH")
  if (is.null(vsearchpath)) vsearchpath <- "~/vsearch-2.22.1/bin/"
  Sys.setenv(PATH = paste(old_path, path.expand(vsearchpath), sep = ":"))
  sample_list <- gsub("_FREYJA_uniq.fasta","",list.files(pattern="^[a-zA-Z0-9]{4}_[a-zA-Z0-9]{4}_sample_[a-zA-Z0-9]{3}_FREYJA_uniq.fasta$"))

  message("HELA will remove chimaeras from each sample")
  X <- NULL
  for (i in sample_list) {
    X <- c(X,paste0("vsearch --uchime_denovo ",i,"_FREYJA_uniq.fasta ",
                    "--sizeout --minh 0.90 ",
                    "--nonchimeras ",i,"_HELA_nonchimeras.fasta ",
                    "--chimeras ",i,"_HELA_chimeras.fasta ",
                    "--uchimeout ",i,"_HELA_uchimeout.log"))
  }
  no_cores <- cores
  clust <- makeCluster(no_cores)
  clusterExport(clust, list("X","old_path","vsearchpath"),envir = environment())
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  stopCluster(clust)

  after_HELA <- mclapply(sample_list,function(file){
    output <- system(paste0("grep '>' ",file,"_HELA_nonchimeras.fasta | wc -l"),intern = T,wait = T)
    value <- as.numeric(output)
    return(data.frame(file=paste0(file,"_HELA_nonchimeras.fasta"),
                      num_seqs=value))
  },mc.cores = cores)
  after_HELA <- do.call("rbind",after_HELA)

  save(file = "summary_HELA.RData",list = c("after_HELA"))

  message("HELA is done.")
}


