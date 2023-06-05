---
title: Alignment for RNA Seq
tags: [formatting]
keywords: Alignment, RNA Seq
last_updated: June 4, 2023 
summary: "Alignment Using HISAT2"
sidebar: mydoc_sidebar
permalink: Alignment_RNA_Seq.html
folder: mydoc
---

# What is Alignment?

Alignment, in the context of sequencing data analysis, refers to the process of matching or aligning short sequence reads obtained 
from high-throughput sequencing technologies to a reference genome or transcriptome. The purpose of alignment is to determine the location 
or position in the reference sequence where each read likely originated from.
In RNA-seq, alignment refers to the process of mapping the RNA sequencing reads obtained from a sample to a reference genome or transcriptome. 
The goal of alignment in RNA-seq is to determine the locations in the reference sequence where the reads likely originated from, enabling the 
identification and quantification of the expressed genes and transcripts.

In this tutorial, alignment is a key step, which provides a significant advantage by enabling gene expression quantification. By aligning the RNA-seq reads to a 
reference genome or transcriptome, each read is assigned to a specific gene or transcript. This allows researchers to accurately quantify the expression levels of 
genes, providing valuable insights into which genes are expressed and at what levels in a given sample. Such information plays a crucial role in understanding 
the underlying biological processes.

# Alignment using HISAT2

## What is HISAT2?

HISAT2 is a popular alignment tool widely used in RNA-seq data analysis. It is specifically designed to align RNA sequencing reads to a reference genome 
or transcriptome accurately and efficiently. HISAT2  enables fast and memory-efficient alignment of reads. HISAT2 employs a hierarchical approach that 
efficiently handles both spliced and unspliced alignments, making it suitable for aligning reads to genomes with a large number of exons and 
alternative splicing events.

A few advantages of HISAT2 are listed below:

  1. Performance and accuracy: HISAT2 is known for its high alignment accuracy and robust performance, even when dealing with large and complex genomes. 

  2. Splice-aware alignment: HISAT2 has the ability to handle spliced alignments, enabling the accurate mapping of reads across exon-exon junctions. 
     This feature is particularly important for studying alternative splicing events and isoform expression analysis.
     
  3. Support for various organisms: Mention that HISAT2 is designed to align reads from a diverse range of organisms, including humans, 
     model organisms such as mouse, and plants. Its compatibility with different species makes it a versatile tool for RNA-seq analysis 
     across various research areas.

## How to use HISAT2.

Here's a simple step-by-step guide on how to use HISAT2 for Alignment:

If you are working on a server, HISAT2 might already be installed. The server may also contains the indexes of the genome of your preference. At the end of the page, you can find instuctions about downloading or building hisat2 indexes de novo(well actually from a fasta genome file).


In order to align the RNA-seq reads to the reference using the HISAT2 aligner:

Use the HISAT2 command with appropriate options to specify the input files, reference index, and other parameters. For example:

For single end data:

      ```
      hisat2 -p [number_of_threads_to_use] -x /path_to_your_hisat2_indexes_folder/genome -U /path_to_your_fastq_files_folder/reads_.fastq -S output.sam
      ```

For paired end data:

      ```
      hisat2 -p [number_of_threads_to_use] -x /path_to_your_hisat2_indexes_folder/genome -1 /path_to_your_fastq_files_folder/reads_R1.fastq -2 /path_to_your_fastq_files_folder/reads_R2.fastq -S output.sam
      ```

When running the HISAT2 alignment, make sure to include the -x hisat2 option, which should point to the directory path and the index_prefix (in this case genome) of your index files (without the .ht2 extension).


   -x This option specifies the path and index prefix of the HISAT2 index files. 
   
   -U This option is used when you have unpaired reads in a single-end sequencing experiment. It specifies the path to the FASTQ file containing the unpaired     reads.
   
   -1 These options are used for paired-end sequencing experiments. -1 specifies the path to the FASTQ file containing the first mate (forward)
   
   -2 These options are used for paired-end sequencing experiments. -2 specifies the path to the FASTQ file containing the second mate (reverse)
   
   -S This option is used to specify the output file where the aligned reads will be written in the SAM (Sequence Alignment/Map) format. 
   
   -p This option is used to specify the number of CPU threads or cores to be used for the alignment process
   
More options about hisat2 can be found [here](http://daehwankimlab.github.io/hisat2/manual/).

The output of hisat2 is a file in SAM format. You can read about this format in [here](https://en.wikipedia.org/wiki/SAM_(file_format)). In brief, this file represents the alignment of sequencing reads to a reference genome or transcriptome. SAM files are typically text-based and can be easily human-readable.      
It's important to refer to the HISAT2 documentation and user manual for detailed information on available options, best practices, and specific considerations for 
your analysis.

## Convert sam files to bam files using samtools

For downstream analysis, we have to convert sam files to bam files. BAM files format is the binary format of SAM files. Converting SAM (Sequence Alignment/Map)
files to BAM (Binary Alignment/Map) format offers several advantages in terms of file size, processing speed, and data storage. This is done by Samtools. Samtools 
is a widely used bioinformatics software package for manipulating and analyzing sequencing alignment files in the SAM (Sequence Alignment/Map) and BAM (Binary 
Alignment/Map) formats. It provides a collection of tools and utilities for tasks such as file format conversion, sorting, indexing, filtering, and statistical 
analysis of alignment data.

      ```
      # Convert sam file to bam file 
      samtools view -b -h -S output.sam > output.bam
      
      # Sort the bam file
      samtools sort output.bam -o output.sorted.bam
      
      # Create an index for the sorted bam file
      samtools index output.sorted.bam
      ```
      
After completing these steps, you can remove the sam file and the unsorted file, in order to free some disk space.


# Further analysis:

Once the alignment is complete and the BAM file is generated, you can proceed with downstream analysis tasks such as gene expression quantification, differential 
expression analysis, variant calling, or other RNA-seq analyses. Various tools and pipelines are available to perform these analyses, such as HTSeq, 
featureCounts, Cufflinks, or DESeq2. Choose the appropriate tool based on your specific research question and analysis goals. In the step, you will learn (hopefully) how to quantify gene expression, using HTSeq. 


# Building indexes for hisat2 

Before aligning reads, you need to build an index of the reference genome or transcriptome using HISAT2. This step needs to be done only once for a given 
reference.

Use the HISAT2-build command to create the index. For example:

```
hisat2-build your_reference_genome.fasta index_prefix
```

Executing this command will generate multiple ht2 format index files with the specified index_prefix.


When running the HISAT2 alignment, make sure to include the -x hisat2 option, which should point to the directory path and the index_prefix of your index files (without the .ht2 extension).

It's worth noting that HISAT2 provides pre-built indexes for various organisms on their [website](http://daehwankimlab.github.io/hisat2/download/). You can conveniently download the index file that matches your data, eliminating the need to build it yourself."




{% include links.html %}
