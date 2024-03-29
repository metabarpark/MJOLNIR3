% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mjolnir2_FREYJA.R
\name{mjolnir2_FREYJA}
\alias{mjolnir2_FREYJA}
\title{FREYJA: Filtering of Reads, Enrollment, Yoke-reads Joining and Alignment}
\usage{
mjolnir2_FREYJA(
  lib_prefix = "",
  cores = 1,
  Lmin = 299,
  Lmax = 320,
  lib = "EXPX",
  fastq_output = T,
  score_obialign = 40,
  demultiplexed = F,
  primer_F = "GGWACWRGWTGRACWNTNTAYCCYCC",
  primer_R = "TANACYTCNGGRTGNCCRAARAAYCA",
  R1_motif = "_R1",
  R2_motif = "_R2",
  obipath = NULL,
  remove_DMS = T
)
}
\arguments{
\item{cores}{Numeric. Number of threads for parallel processing.}

\item{Lmin}{Numeric. Minimum bp length for a sequence to be accepted.}

\item{Lmax}{Numeric. Maximum bp length for a sequence to be accepted.}

\item{lib}{Character string. Acronym for the experiment. This
acronym must be of 4 characters in capital letters. Do not mix up library and
experiment acronyms. However they can be the same.}

\item{demultiplexed}{Logical. If TRUE, the data has been already demultiplexed. 
otherwise this has to be set to FALSE. See details for further information.}

\item{primer_F}{Character string of the Forward primer. Necessary when 
demultiplexed=TRUE.}

\item{primer_R}{Character string of the Reverse primer. Necessary when 
demultiplexed=TRUE.}

\item{R1_motif}{Character string that distinguish the forward line file from
the reverse.}

\item{R2_motif}{Character string that distinguish the reverse line file from
the forward.}

\item{obipath}{Character string specifying the PATH to the obi binary.}

\item{remove_DMS}{Logical. If TRUE, it will delete all obidms objects that are
created during the process. This can save a lot of hard disk space. The FALSE 
option is useful for developing and debugging.}

\item{lib_prefixes}{Character vector. Acronym for each sequencing library. This
acronym must be of 4 characters in capital letters. Do not mix up library and
experiment acronyms. The latter will be required in following steps. However 
they can be the same.}

\item{fastq_ouput}{Logical. If TRUE, a fastq file for each sample with the 
demultiplexed and quality filtered sequences will be retrieved. This is useful
 if you wish to publish your data in public repositories.}

\item{score_obilign}{Numeric. Minimum quality threshold to retain a sequence 
in the qualiy filtering after pairalignment.}
}
\description{
FREYJA will use OBITools3 commands to merge paired-end reads, demultiplex 
libraries into samples (if needed),trim primer sequences, filter by length, 
split sequences per sample and dereplicate within each sample.
}
\details{
In case the data is already demultiplexed and consist of individual fastq 
files for each sample, use the option demultiplexed=TRUE. When demultiplexed=TRUE, 
FREYJA will read the names of each individual R1 fastq files from a column 
in the LIBX_metadata.tsv file, called fastq_name_R1.
In the metadata table, each sample in the original_samples column must have 
a matching fastq_name_R1 and a matching mjolnir_agnomen (LIBX_sample_XXX).
When demultiplexed=TRUE, you must also specify the primer_F and primer_R 
sequences in the options input to FREYJA. COI Leray-XT primers are specified 
by default. Otherwise, when demultiplexed=FALSE, the primers information 
must be already written in the LIBX_ngsfilter.tsv files.
Until further optimization, please do not set more than 7 cores for the parallel
processing of the next step, FREYJA.
}
\examples{
library(mjolnir)

# Define input fastq files (only names of R1 files are needed)
R1_filenames <-c("ULO1_R1.fastq.gz","ULO2_R1.fastq.gz","ULO3_R1.fastq.gz","ULO4_R1.fastq.gz")

# Input identifiers for the individual libraries to be used. It should be a 4-character name, matching the information in the ngsfilter files
lib_prefixes <- c("ULO1","ULO2","ULO3","ULO4")

# Enter number of cores to be used in parallel. 
cores <- 7

# Run RAN
mjolnir1_RAN(R1_filenames,cores,lib_prefixes,R1_motif="_R1",R2_motif="_R2")

# set experiment acronym
lib <- "ULOY"

# Run FREYJA
mjolnir2_FREYJA(lib_prefix = lib_prefixes,lib = lib,cores = cores,Lmin=299,Lmax=320)
}
