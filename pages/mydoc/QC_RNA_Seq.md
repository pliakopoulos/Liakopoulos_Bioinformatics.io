---
title: Tooltips
tags: [formatting]
keywords: popovers, tooltips, user interface text, glossaries, definitions
last_updated: June 4, 2023 
summary: "Quality Control Using FastQC"
sidebar: mydoc_sidebar
permalink: Create QC_RNA_Seq.html
folder: mydoc
---

## What is Quality contol?

Quality control (QC) involves assessing the quality, integrity and evaluation of the sequencing data obtained from RNA-seq experiments. It involves various 
steps to identify and address potential issues that can affect the accuracy and reliability of downstream analysis.The goal is to identify any potential issues 
that could affect downstream analysis and to ensure that the data is reliable and of high quality. 

The key aspects of pre - alignment quality control for RNA - seq are:

  1. Sequence Quality Assessment: This step involves evaluating the quality of individual sequencing reads. Quality scores, represented as Phred scores, 
     indicate the confidence level of each base call. QC tools such as FastQC or Fastp can be used to assess parameters such as read length distribution, 
     per-base sequence quality, GC content, and the presence of sequencing adapters.

  2. Adapter and quality trimming: If adapter sequences are present in the reads, or if there are low-quality regions, it is important to trim them. 
     Trimming removes adapter sequences and low-quality bases, improving the overall quality of the data and reducing the influence of sequencing artifacts.

## Quality Control Using FastQC

# What is FastQC?

FastQC is a widely used tool for performing quality control (QC) analysis on raw sequencing reads obtained from RNA-seq experiments. As explained in their [documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/1%20Introduction/1.1%20What%20is%20FastQC.html),
FastQC aims to provide a QC report which can spot problems which originate either in the sequencer or in the starting library material.

# How to use FastQC

Here's a simple step-by-step guide on using FastQC for QC analysis:

  1. Prepare input files: Make sure you have the raw sequencing data in a format supported by FastQC, such as FASTQ or BAM. It also supports GZip compressed FastQs.

  2. Run FastQC: To run FastQC in your terminal you only need to run:
      {% raw %}
      ```
      ## Run FASTQC
      fastqc -o [path_to_your_output_directory] -t [number_of_threads_to_use] input_file.fastq
      ```
      {% endraw %}
     The -o and -t flags refer to:
     
      -o Create all output files in the specified output directory.
     
      -t Specifies the number of files which can be processedsimultaneously.  
         Each thread will be allocated 250MB of memory so you shouldn't run more threads than your available memory will cope with.
      
     Using FastQC with -h flag and withoun any other option, the help file of fastqc will be printed. This file contains information about the different flags 
     this tool possesses. It also includes information about how to correctly run the programm, and a sort description.
   
   3. Evaluation of the results.

More information about FastQC can be found it this manual [website](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf)

# Evaluating Results
   
   
   

Interpret the results: After the analysis is complete, FastQC generates a comprehensive report that provides various metrics and visualizations 
for assessing the quality of the sequencing data. 
Open the HTML report(s) to view the results.

# FastQC outputs

The FastQC report contains several sections, including:

Basic Statistics: This section provides general information about the dataset, such as the total number of reads, sequence length distribution, and GC content.

Per-base Sequence Quality: Here, FastQC displays the quality scores for each base position in the reads. 
It helps identify regions with low-quality scores, potential sequencing errors, or biases.

Per-sequence Quality Scores: This section presents a graph showing the distribution of average quality scores for all reads. 
It allows you to assess the overall quality of the dataset.

Sequence Length Distribution: This plot shows the distribution of read lengths in the dataset, helping identify any biases or unexpected read truncations.

Overrepresented Sequences: FastQC identifies sequences that occur at an unusually high frequency, which could indicate contaminants or adapter sequences that 
need to be trimmed.

Adapter Content: This section helps identify if adapter sequences are present in the reads, which may need to be removed during preprocessing.

Interpret the warnings and flags: FastQC flags potential issues and generates warnings if certain metrics fall outside the expected range. 
Pay attention to these warnings as they highlight potential problems with the data, such as low-quality scores, overrepresented sequences, 
or adapter contamination.
By examining the FastQC report, you can gain insights into the quality of your sequencing data 
and take appropriate measures for preprocessing or downstream analysis.





## Adapter trimming and low reads filtering using TrimmGalor








{% raw %}
```html
<a href="#" data-toggle="tooltip" data-original-title="{{site.data.glossary.jekyll_platform}}">Jekyll</a> is my favorite tool for building websites.
```
{% endraw %}

This renders to the following:














{% include links.html %}
