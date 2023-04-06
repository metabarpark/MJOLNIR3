# MJOLNIR3
<p align="center">
  <img src="https://github.com/adriantich/MJOLNIR3/blob/main/MJOLNIR.png">
</p>

<H1><b>Metabarcoding Joining Obitools &amp; Linkage Networks In R</b></H1>

<b>by Owen S. Wangensteen & Adrià Antich.</b>

MJOLNIR3 is a powerful tool to crush big amounts of raw metabarcoding data, and molding them into organized data sets of taxonomically assigned MOTUs. 

MJOLNIR3 comes in an R package, so that modular metabarcoding pipelines are easy to run from the R environment. MJOLNIR3 runs on Linux and Mac systems. The extensive use of package parallel and several dependencies that are designed primarily for Linux systems (see below) makes the success of installations in Windows highly improbable. Users are welcome to try to install and run MJOLNIR3 on Windows Linux Subsystem, but we would not recommend that.

MJOLNIR3 depends on the following dependencies, which must be installed in the system and properly working:

- OBITools3:\
  This is the new version of the Obitools 2 ([Boyer et al. 2016](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12428)) \
  Original information about OBITools here: https://metabarcoding.org/obitools3 \
  Help on installing OBITools: http://rleca.pbworks.com/w/file/fetch/124098201/tuto_obitools_install_W10OS.html \

- VSEARCH ([Rognes et al. 2016](https://peerj.com/articles/2584/)) \
  Help on installing VSEARCH: 
  https://github.com/torognes/vsearch
  
- DnoisE ([Antich et al. 2022](https://peerj.com/articles/12758/)) \
  Help on installing DnoisE: 
  https://github.com/adriantich/DnoisE
  
- SWARM v2.0 or newer ([Mahé et al. 2015](https://peerj.com/articles/1420/)) \
  Help on installing SWARM: 
  https://github.com/torognes/swarm
  
- LULU ([Frøslev et al. 2017]()) \
  Help on installing LULU:
  https://github.com/tobiasgf/lulu

- Package Biostrings from the Bioconductor suite. \
  Help on installing Biostrings:
  https://bioconductor.org/packages/release/bioc/html/Biostrings.html
                   

### Installing MJOLNIR3:

1. Create and activate a virtual environments: MJOLNIR3 is highly recommended to be run in python environments. Python3.6 or higher is required but the user can choose its best option. Options: 
  
- [venv](https://docs.python.org/3/library/venv.html) (sudo required to install it) 

- [pyen](https://github.com/pyenv/pyenv)

- conda ([Anaconda](https://docs.anaconda.com/anaconda/install/index.html) or [miniconda](https://docs.conda.io/en/latest/miniconda.html))

2. Install dependencies with the virtual environment activated. With the Conda environment some of the software can be installed from conda repositories.

* It is recommended to download the repository, so constant updates will be occurring for the first versions and updates are easier that way. Then install all the required software within the MJOLNIR3 folder.

        # with your prefered environment activated
        # clone MJOLNIR3 repository
        git clone https://github.com/metabarpark/MJOLNIR3.git
        cd MJOLNIR3
        # installation of Obitools3
        git clone https://git.metabarcoding.org/obitools/obitools3.git
        cd obitools3
        pip3 install cython
        python3 setup.py install
        source obi_completion_script.bash
        cd ..
        # installation of vsearch
        git clone https://github.com/torognes/vsearch.git
        cd vsearch
        ./autogen.sh
        ./configure CFLAGS="-O3" CXXFLAGS="-O3"
        make
        cd ..
        # installation of swarm
        git clone https://github.com/torognes/swarm.git
        cd swarm/
        make
        cd ..
        # installation of DnoisE
        git clone https://github.com/adriantich/DnoisE.git
        cd DnoisE/
        python3 setup.py install
        cd ../..
        # now turn to R
        R
        # install lulu
        > library(devtools)
        > install_github("tobiasgf/lulu")  
        # install biostrings
        > if (!require("BiocManager", quietly = TRUE))
               install.packages("BiocManager")
        > BiocManager::install("Biostrings")
        # Finally install MJOLNIR3
        > install.packages('MJOLNIR3', repos = NULL)

* For any update remove the package and install it again. In R console:

        > remove.package('mjolnir')
        > install.packages('MJOLNIR3', repos = NULL)

* Alternatively you can install each dependency using anaconda or your prefered method but keep in mind that the paths to these dependencies will have to be specified in the conda paths, in your own PATH or specified for each in your script.

3. Alternative ways to install MJOLNIR3 
* If the devtools library is properly installed in the system: MJOLNIR3 can be installed directly from the R console using: 

        > library(devtools)
        > install_github("metabarpark/MJOLNIR3") 
* If the devtools library is not installed: 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; MJOLNIR3 can be downloaded as a package from this link: https://github.com/metabarpark/MJOLNIR3/archive/main.zip.
Then the file must be unzipped and MJOLNIR3 can be installed offline from the R console using:

        > install.packages("MJOLNIR3-main", repos=NULL)
        
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Another option is to clone with git the MJOLNIR3 repository: 

        git clone https://github.com/metabarpark/MJOLNIR3.git
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; and then install the package from R:

        > install.packages("MJOLNIR3-main", repos=NULL)

MJOLNIR3 is currently optimized by default to process COI metabarcoding data (the Leray Fragment). For other metabarcoding markers, some settings must be changed as follows:

The following settings are recommended for COI Leray/Leray-XT primers (Leray et al. 2013; Wangensteen et al. 2018): 
- In mjolnir2_FREYJA: Lmin=299,Lmax=320 
- In mjolnir4_ODIN: d=13

The following settings are recommended for 12S MiFish primers (Miya et al. 2015): 
- In mjolnir2_FREYJA: Lmin=140,Lmax=190 
- In mjolnir4_ODIN: d=1,algorithm="SWARM"

The following settings are recommended for 18S All_shorts primers (Guardiola et al. 2015): 
- In mjolnir2_FREYJA: Lmin=75,Lmax=180 
- In mjolnir4_ODIN: d=1,algorithm="SWARM"

The following settings are recommended for 16S Bacterial F515/R806 primers (Caporaso et al. 2011): 
- In mjolnir2_FREYJA: Lmin=215,Lmax=299 
- In mjolnir4_ODIN: d=1,algorithm="SWARM"

<H2>The MJOLNIR3 Pipeline</H2>

This is a simplified scheme of the MJOLNIR3 workflow:

<p align="center" width="400">
  <img src="https://github.com/adriantich/MJOLNIR3/blob/main/workflow_main.png">
</p>

<B>0. Input data</B>

MJOLNIR3 can process paired-end FASTQ files from Illumina sequencing. There are two ways that you can use to start MJOLNIR3: 

<B>[1] The MAGENTA WAY</B> (in the figure) will be used if you have <B>MULTIPLEXED</B> libraries containing information from several samples (typically prepared using the METAFAST protocol). This protocol adds sample-tags on both ends of the amplicons, at 5' from the metabarcoding primers. Many samples can be multiplexed into a single metafast library. Each pair of fastq files belonging to a library usually contain multiple samples, identified by unique combination of forward:reverse sample-tags. MJOLNIR3 can process several such libraries simultaneously, spanning hundreds or thousands of samples, which will be joined together into a final combined dataset before the clustering step. Then you need to start your pipeline in mjolnir1_RAN(), so you can increase the speed using parallel sequencing. You will need to provide ngsfilter files for each library, that will be used by mjolnir2_FREYJA(). You have an example of a pipeline using this way here: https://github.com/metabarpark/MJOLNIR3/tree/main/example_MJOLNIR_multiplexed_data_metafast

<B>[2] The ORANGE WAY</B> (in the figure) will be used if you have individual <B>DEMULTIPLEXED</B> paired-end fastq files for each sample. That is, if the demultiplexing have been already done by the Miseq analyzer (typically using library tags). These fastq files are the raw files from the sequencer, and must not be pre-trimmed, so they still contain the primer sequences. Then you will go directly to the second step, using the option demultiplexed=TRUE: mjolnir2_FREYJA(demultiplexed=T). You will need to provide the list of fastq files to process in the EXPX_metadata.tsv table. You have an example of a pipeline using this way here: https://github.com/metabarpark/MJOLNIR3/tree/main/example_MJOLNIR_demultiplexed_data

<B>1. RAN: Reads Allotment in N portions </B>

MJOLNIR3 starts with a call to mjolnir1_RAN() which will split initial FASTQ files in aliquot parts for parallel processing. This will distribute the initial paired-end fastq files in fragments of equal number of reads, to be processed in parallel. You do not need to worry about if your fastq files are gzipped (.fastq.gz) or unzipped (.fastq). RAN will be able to process both formats automatically. You have to specify the R1 and R2 motifs, needed to get the name of the R2 fastq file from its R1 counterpart. (WARNING! in this step it is recommended to set cores in parallel mode to no more than 7 for a better hard disk memory management in FREYJA. This will be optimized soon.)

<B>2. FREYJA: Filtering of Reads, Enrollment, Yoke-reads Joining and Alignment </B>

The function mjolnir2_FREYJA() calls OBITools3 for the following three steps: (1) paired-end alignment, (2) demultiplexing and (3) read-length filtering, wrapped together in a single function. Paired-end alignment is done using alignpairedend. Removal of low-quality reads is based on the individual read quality results from the latter step. Demultiplexing and primer-removal are done using ngsfilter. A length filter using grep follows, based on the length of the metabarcoding fragment (excluding the primers). All reads having bases different from [ACGT] are also removed.
In case you are using the magenta way, the ngsfilter function needs an input table (ngsfilter-table), containing information about the internal sample names (agnomens), the sample-tags used to identify them at both ends, and the sequences of the metabarcoding primers used. ngsfilter-tables must be provided for each library. You have some examples of ngsfilter files in the <A href="https://github.com/github.com/metabarpark/MJOLNIR3/tree/main/Example_Pipeline_Metafast">Example Pipeline Metafast folder</A>. Note that the internal names of the samples (agnomens) written in the ngsfilter files MUST have a particular format, beginning with a 4-character library code (the same library code must be specified within the lib_prefixes variable at the beginning of the pipeline), and ending in a three-digit numerical code; e.g.: **LIBX_sample_XXX**. Also, the name of the ngsfilter file for each library has to follow a fixed syntax: **ngsfilter_LIBX.tsv**. They must be tab-separated text files, following the same format than is used in OBITools command ngsfilter (https://pythonhosted.org/OBITools/scripts/ngsfilter.html). 
In case you are using the orange way, you have to use the option demultiplexed=TRUE. Then you do not have to provide any ngsfilter table. But you must specify the names of the individual fastq files for each sample in the EXPX_metadata.tsv file, under a column called "fastq_name_R1". In this case, you have to specify the R1 and R2 motifs here (since you are not using mjolnir1_RAN), and you have to provide the sequences of the forward and reverse primers (the COI Leray-XT sequences are provided by default).
In both cases, the output of FREYJA will be individual, aligned, filtered, trimmed, fasta files, containing exclusively the sequence of the metabarcoding marker, with the format EXPX_LIBX_sample_XXX_FREYJA.fastq with the option fastq_ouput=TRUE. One such file for each sample. Dereplication of sequences is then performed with the uniq function from Obitools3 within each sample retrieving a fasta file for each sample with the name EXPX_LIBX_sample_XXX_FREYJA_uniq.fastq
(WARNING! in this step it is recommended to set cores in parallel mode to no more than 7 for a better hard disk memory management. This will be optimized soon.)
**LIBX refers to each library acronym (lib_prefix parameter)** & **EXPX refers to the experiment acronym (lib parameter)** and they can be the same

<B>3. HELA: Hierarchical Elimination of Lurking Artifacts</B>

MJOLNIR3 uses the mjolnir3_HELA() function to call the uchime_denovo algorithm implemented in VSEARCH, to remove chimaeric sequences from the individual sample files provided by FREYJA, in a sample-by-sample basis. This procedure is much faster than removing the chimaeras from the whole dataset together. This procedure ends with a fasta file for each sample with the name EXPX_LIBX_sample_XXX_HELA_nonchimeras.fasta. This step can be parallelized using the maximum cores available in the computing machine but no more that the number of samples itself.

<B>4. ODIN: OTU Delimitation Inferred by Networks</B>

The function mjolnir4_ODIN() uses the four different strategies to delimit MOTUs and/or ESVs. This strategies are set with the algorithm parameter: a)"DnoisE_SWARM", b)"SWARM", c)"SWARM_DnoisE" and d)"DnoisE". In short, DnoisE refers to the denoising process with DnoisE to obtain ESV and SWARM to a clustering process with SWARM to obtain MOTUs. 
DnoisE is a software to merge spurious sequences into their "mothers" (see Antich et al. 2022) to obtain Exact (also Amplicon) Sequence variants. DnoisE is an open source and parallelizable alternative to Unoise that allows to introduce an entropy correction based on the different entropies of each position in the codon of coding genes. This is highly recommended for markers as COI for which this program was intended. However, with the entropy=FALSE parameter, this programs performs the same denoising procedure as described for Unoise3.
SWARM is an algorithm to delimit MOTUs, based on linkage-networks created by step-by-step agregation. This clustering algorithm is not based on an arbitrary, constant, absolute identity threshold. Conversely, SWARM is based on an iterative aggregation of sequences that differ less than a given distance d. This strategy results into linkage-networks of different sizes, which can have different effective values for within-MOTU identity threshold, depending on the complexity of the natural variability of the sequences present in the sample. This procedure is very convenient in the case of hypervariable metabarcoding markers, such as COI, which usually feature extremely high levels of natural diversity, in addition to the random sequencing errors. 

## WARNING!! From here and below, the README.md file is not updated. Please do not take this into account, we are working on it!

After removing the chimaeras from every sample, HELA joins all the samples again, and de-replicates the sequences. Then creates new sequence identifiers for all different sequence variants found, and generates two outputs: (1) a fasta file with total abundances (ending in .vsearch.fasta) ready to be fed into Swarm for the clustering step, and (2) a table including the abundance information of all sequence variants in every sample. 

Moreover, the resulting networks are very convenient for further processing of metaphylogeography datasets (Turon et al. 2019). In the future, mjolnir4_ODIN() will be able to produce ESV tables (a proxy for haplotypes) using the option generate_ESV=TRUE and several option to denoising the data. However, this option is still in beta test and it has not been properly implemented yet. 
mjolnir4_ODIN() takes the LIBX.vsearch.fasta file produced by HELA, as an input. Then it uses the table with read abundances by sample to compute the combined abundance of the resulting MOTUs in every individual sample. The output is a CSV file with contains a row for each MOTU, with its abundances in every sample. Also, a fasta file including the representative sequences of each MOTU (usually only non-singleton MOTUs are included) is generated, which will be used for taxonomic assignment in the next step.

<B>5. THOR: Taxonomy with Higher-than-Order Ranks</B>

Taxonomic assignment is performed with the mjolnir5_THOR() function, which is a wrapper of ecotag (Boyer et al. 2016) and owi_add_taxonomy (Wangensteen & Turon 2017). This step depends on the availability of a taxonomic database in ecoPCR format (from which the phylogenetic tree information is retrieved) and a reference database file including sequences for the exclusively metabarcoded fragment, conveniently identified by a taxonomic identifier (taxid), in fasta format. These can be downloaded from our DUFA repository: https://github.com/uit-metabarcoding/DUFA. mjolnir5_THOR then uses a custom procedure to assign higher taxonomic ranks (higher than order). These taxonomic ranks are stored in a csv file called "order_complete.csv", which is customizable to meet the particular preferences of the user. This is specially useful for assigning the higher ranks of microeukaryotes, whis is remarkably unstable and there is not universally-accepted consensus. 

<B>6. FRIGGA: Final Recount and Integration of Generated Genealogies and Abundances</B>

The function mjolnir6_FRIGGA() will join the taxonomic information obtained by THOR with the information of abundances per sample calculated by ODIN. FRIGGA will create an output CSV file ending in .All_MOTUs.csv. In the future, FRIGGA will have a more important role in integrating MOTU-level data with ESV-level data (species and haplotypes), but this option still needs implementation. 

<B>7. LOKI: LULU Overseeing with Kinship Identification </B>

The next step of MJOLNIR3 is the removal of pseudogenes using the mjolnir7_LOKI() function. This is a wrapper of the LULU algorithm. The information about the removed putative pseudogene MOTUs is also stored as an additional output file, together with the information of their possible mother sequences. This file can be checked to assess the taxonomic coherence of the results. A .curated.csv file is created.

<B>8. RAGNAROC: Replace AGnomens with Names And Recover Original Codification</B>

mjolnir8_RAGNAROC() is the last step of the MJOLNIR3 pipeline. It requires a tsv metadata table (default name: LIBX_metadata.tsv) with a column of the original sample names and another column with the internal agnomens used along the rest of the pipeline. RAGNAROC also reorders the sample columns by alphabetical order and produces the LIBX.final_dataset.csv file. RAGNAROC will apply a Relative Read Abundance (RRA) filter on a sample by sample basis, removing reads of MOTUs that are below a relative RRA threshold (default min_relative=1/50000) for every particular sample. RAGNAROC will then apply a total abundance filter to the MOTUs, removing MOTUs that have lower total abundance than a given value (default min_reads=2). RAGNAROC will also optionally remove those MOTUs that have been assigned to prokaryotic sequences or to the root of the Tree Of Life. Those MOTUs that are considered contaminants (as specified in a dedicated text file) can also be automatically removed. The MOTUs surviving RAGNAROC can also be ordered by their taxonomy or their abundance before writing the LIBX.final_dataset.csv file.
