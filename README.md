# circRNAFull: an R package for reconstruction of full length circRNA sequence using chimeric alignment information

## Requirements
* R (>= 4.1.2)
* Biostrings
* seqinr
* stringi
* Rsamtools

## Installation
### From cran
To install the package from cran, run the command:

    install.packages("circRNAFull", dep=T)
	
### From github
To install the package from github first you need to install the package “devtools” using the following command:

    install.packages("devtools", dep=T)

The package "circRNAFull" depends on two bioconductor packages "Biostrings" and "Rsamtools" which cannot be installed automatically while installing "circRNAFull" using "devtools". So, you need to install these two packages manually using the following way:

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("Biostrings")
    BiocManager::install("Rsamtools")

Finally, install “circRNAFull” by the following command:

    devtools::install_github("tofazzal4720/circRNAFull", dep = T)

Start analysis by loading the package with the following command:

    library("circRNAFull")

## Extracting transcript name and spanning reads for circRNAs

#### Description
This function makes a text file containing circle ids, transcript name and the spanning redas for each circRNAs.

#### Usage
    
    extract_circle_ids(circ_out, chimeric_out, output)
    
#### Arguments

`circ_out` is the output from the circRNA prediction tool CIRCexplorer.

`chimeric_out` is the junction file obtained from the chimeric alignment produced by STAR aligner.

`output` is the name of the output folder.

#### Value
The circle ids, transcript name and the spanning reads will be produced in '*output*' folder.

#### Example
    
    #Loading an example output from the circRNA prediction tool CIRCexplorer
    data(CIRCexplorer_output)
    circ_out<-CIRCexplorer_output
    
    #Loading an example junction file obtained from the. chimeric alignment produced by STAR
    data(Chimeric.out.junction)
    chimeric_out<-Chimeric.out.junction
    
    #A temporary directory is created as an output folder.
    output<-tempdir()
    
    #Extracting circle ids, transcript name and spanning reads of circRNAs. 
    #The output will be written in 'output' directory
    extract_circle_ids(circ_out, chimeric_out, output)

## Extracting individual alignment file for each circRNA from chimeric alignment bam file
#### Description
This function generates individual alignment file for each circRNA from the chimeric alignment produced by STAR.

#### Usage
    extract_reads_from_bam<-function(circle_id, bamfile, outfolder)

#### Arguments
`circle_id` is a data frame containg circle ids, transcript name and spanning reads of circRNAs obtained from function
`extract_circle_ids()`.

`bamfile` is a chimeric alignment bam file produced by STAR read as a data frame.

`outfolder` is the name of the output directory.

#### Value 
Individual alignment files for each circRNAs will be produced in *outfolder*.

#### Example
    #loading an example circle_id file
    data(circle_ID)
    circle_id<-circle_ID
    
    #Please upload your chimeric alignment bam file.Suppose you have the chimeric
    #alignemnt bam file 'Chimeric.out.bam' in you working directory.
    #Then run:
    bam <- scanBam("Chimeric.out.bam")
    bamfile <- as.data.frame(bam)
    
    #Creating an output directory
    outfolder<-tempdir()
    
    #Extracting individual alignment file for each circRNAs. 
    #The individual alignment file will be generated in the 'outfolder' directory.
    extract_reads_from_bam<-function(circle_id, bamfile, outfolder)
 
 ## Extracting exon for the circRNAs
 
 #### Description
 This function extracts exons for the circRNAs from the CIGAR value of the chimeric alignment obtained from *STAR*.
 
#### Usage
    Extract_exon(bed_file, folder_name, circle_ids, output_folder)
    
 #### Arguments
 
 `bed_file` is a data frame containing the exon coordinates obtained from reference genome annotation.
 
 `folder_name` is the location of the individual alignments obtained from function `extract_reads_from_bam()`.
 
 `circle_ids` is a data frame containing the circle id, transcript name and spanning reads obtained from function `extract_circle_ids()`.
 
 `output_folder` is the output folder name.
 
 #### Value
 A file *exon_count.txt* containing exons for each circRNA will be produced in *output_folder*.
 
 #### Examples
 
    #Loading an example bed_file
    data(exon_coordinate)
    bed_file<-exon_coordinate
    
    #Example of folder_name containing the individual alignments
    folder_name<-system.file('extdata', package = 'circRNAFull')
    
    #Loading an example circle_ids file
    data(circle_ID)
    circle_ids<-circle_ID
    
    #creating a temporary output_folder
    output_folder<-tempdir()
    
    #Extracting exons for all circRNAs
    Extract_exon(bed_file, folder_name, circle_ids, output_folder) 
    
## Detecting skip exon
#### Description

This function detects skip exon.

#### Usage

    skip_exon_detection(circle_reads_folder, exon_count_file, output_folder)

#### Arguments

`circle_reads_folder` is the location of the individual alignments obtained from function `extract_reads_from_bam()`.

`exon_count_file` is a data frame containing exons for circRNAs obtained from function `Extract_exon()`.

`output_folder` is the name of the output directory.

#### Value
 A file *detect_skip_exon.txt* containing the skip exons will be produced in *output_folder*.
 
#### Examples
    #Example of circle_reads_folder containing the individual alignments
    circle_reads_folder<-system.file('extdata', package = 'circRNAFull')
    
    #Loading an example exon_count_file containing exons for circRNAs
    data(exon_count)
    exon_count_file<-exon_count
    
    #Creating a temporary output directory
    output_folder<-tempdir()
    
    #Detecting skip exon
    skip_exon_detection(circle_reads_folder, exon_count_file, output_folder)
   
## Extracting exon after deleting skip exon
#### Description
This function extracts the remaining exons after deleting skip exons.

#### Usage
    extract_exon_after_delete_skip_exon(exon_count_file, skip_exon_file, output_folder)

#### Arguments

`exon_count_file` is a data frame containing exons for circRNAs obtained from function `Extract_exon()`.

`skip_exon_file` is a data frame containing skip exons for circRNAs obtained from function `skip_exon_detection()`.

`output_folder` is the name of the output directory.

#### Value
A text file *extract_exon_after_skipping_exon.txt* containing the exons after deleting skip exon will be produced in *output_folder*.

#### Examples
    #Loading an example exon count file
    data(exon_count)
    exon_count_file<-exon_count
    
    #Loading an example skip exon file
    data(skip_exon)
    skip_exon_file<-skip_exon
    
    #Creating a temporary output directory
    output_folder<-tempdir()
    
    #Extracting exon after deleting skip exon
    extract_exon_after_delete_skip_exon(exon_count_file, skip_exon_file, output_folder)
    
## Reconstruction of full length circRNA sequences
#### Description
This function reconstructs the full length circRNA sequences.

#### Usage
    extract_sequence(seq_name, sequence, exon_count_file, output)
#### Arguments
`seq_name` is a vector containing the name of chromosomes of the reference genome.

`sequence` is a vector containing the sequences of the reference genome.

`exon_count_file` is the exon count file after deleting skip exons.

`output` is the output folder name.

#### Value
The reconstructed circRNA sequences *circRNAseq.fasta* will be produced in folder *output*.

#### Examples
    #Please download the reference genome. 
    #Suppose you have downloaded the reference genome 'hg38.fa' in you working directory.
    #Then run:
    ref_genome<-"hg38.fa"
    fastaFile <- readDNAStringSet(ref_genome)
    seq_name = sub('\\ .*', '', names(fastaFile))
    sequence = paste(fastaFile)
    
    #Loading an exon count file after deleting skip exon
    data(exon_count_after_skipping_exon)
    exon_count_file<-exon_count_after_skipping_exon
    
    #Creating an output directory
    output<-tempdir()
    
    #Reconstructing the circRNA sequences. 
    #The circRNA sequences will be generated in the 'output' directory.
    extract_sequence(seq_name, sequence, exon_count, output)
    
