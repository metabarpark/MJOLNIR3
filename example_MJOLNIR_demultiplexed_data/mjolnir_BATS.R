# Load MJOLNIR silently
suppressPackageStartupMessages(library(mjolnir))

# Define number of cores to be used in parallel.
cores <- 16

obipath="~/obi3-env/bin/"
vseatchpath = "~/vsearch-2.22.1/bin/"
swarmpath = "~/swarm/bin/"
dnoise_path = "~/DnoisE/src/"
tax_dir = "~/taxo_NCBI/"
tax_dms_name = "DUFA_COI"

# Input name for the final combined library (should be a 4-character name)
lib <- "BATS"

####################
# MJOLNIR pipeline #
####################


# FREYJA will do the paired-end alignment, demultiplexing & length filtering. This will be done for each sample file separately.
mjolnir2_FREYJA("", cores, Lmin=299, Lmax=320,lib,demultiplexed=T,primer_F="GGWACWRGWTGRACWNTNTAYCCYCC",primer_R="TANACYTCNGGRTGNCCRAARAAYCA",
                R1_motif="_R1.",R2_motif="_R2.",fasta_output = T,fastq_output = T,remove_DMS = T,obipath = obipath)

# HELA will remove chimaeric sequences in a sample-by-sample basis, will change identifiers of remaining unique sequences
# And will generate a table of their abundances in each sample & a fasta file with unique sequences and their total abundance for ODIN
mjolnir3_HELA(lib, cores, vsearchpath = vsearchpath)

# ODIN will do the clustering & will generate a table with the abundances of each MOTU in each sample
mjolnir4_ODIN(lib, cores,d=13,min_reads_MOTU=2,min_reads_ESV=2,alpha=5,COI=T,entropy=c(0.47,0.23,1.02,313),algorithm="DnoisE_SWARM",obipath=obipath, swarmpath=swarmpath, dnoise_path=dnoise_path, remove_singletons = TRUE,remove_DMS=T)

# THOR will asign the taxonomy to the representative sequence of each MOTU
mjolnir5_THOR(lib, cores, tax_dir=tax_dir, tax_dms_name=tax_dms_name, obipath=obipath, run_ecotag=T,remove_DMS=T)

# FRIGGA will integrate the information of MOTU abundances and taxonomy assignment from ODIN & THOR in a single table
mjolnir6_FRIGGA(lib)

# LOKI kill remove the pseudogenes and will keep track of the taxonomic information of the removed MOTUs
mjolnir7_LOKI(lib, min_id=.84,vsearchpath = vsearchpath)

# RAGNAROC will change the names of the samples to recover the original names and will remove unnecessary columns
mjolnir8_RAGNAROC(lib,min_reads=2,min_relative=1/50000,
		remove_bacteria=T,remove_contamination=F,
                ESV_within_MOTU=T,blank_col="BLANK",blank_tag=T,remove_numts=F,cores=1,blank_relative=0.1)
