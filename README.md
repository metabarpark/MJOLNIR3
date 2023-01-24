# MJOLNIR 3
![MJOLNIR](https://github.com/uit-metabarcoding/MJOLNIR/blob/main/mjolnir_blue_mid.png)

<H1><b>Metabarcoding Joining Obitools &amp; Linkage Networks In R</b></H1>

<b>by Owen S. Wangensteen (University of Barcelona).</b>
<b>and Adrià Antich (CEAB-CSIC).</b>

MJOLNIR is a powerful tool to crush big amounts of raw metabarcoding data, and molding them into organized data sets of taxonomically assigned MOTUs. 

MJOLNIR comes in an R package, so that modular metabarcoding pipelines are easy to run from the R environment. MJOLNIR runs on Linux and Mac systems. The extensive use of package parallel and several dependencies that are designed primarily for Linux systems (see below) makes the success of installations in Windows highly improbable. Users are welcome to try to install and run MJOLNIR on Windows Linux Subsystem, but I would not recommend that.

MJOLNIR depends on the following dependencies, which must be installed in the system and properly working:

- OBITools 3 (Boyer et al. 2016):
  Original information about OBITools 3 here: https://git.metabarcoding.org/obitools/obitools3

- VSEARCH (Rognes et al. 2016): 
  Help on installing VSEARCH: https://github.com/torognes/vsearch
  
- SWARM v2.0 or newer (Mahé et al. 2015):
  Help on installing SWARM: https://github.com/torognes/swarm
  
- LULU (Frøslev et al. 2017):
  Help on installing LULU:
  https://github.com/tobiasgf/lulu

- Package Biostrings from the Bioconductor suite. https://bioconductor.org/packages/release/bioc/html/Biostrings.html 
  To install Biostrings, type the following commands in the R console:
 
      if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")
      BiocManager::install("Biostrings")
