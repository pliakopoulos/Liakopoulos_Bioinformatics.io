---
title: About Ruby, Gems, Bundler, and other prerequisites
tags: [RNA Seq, Introduction]
keywords:
summary: "RNA Seq is an expirement for the identification of the gene expression profiles of cells."
sidebar: mydoc_sidebar
permalink: RNA_Seq_Introduction.html
folder: mydoc
---

## RNA Seq

RNA-Seq (RNA sequencing) is a high-throughput technique used to study gene expression at the transcriptome level. It involves the sequencing and analysis of RNA
molecules present in a sample, providing valuable information about the types and quantities of RNA transcripts produced by an organism or tissue at a given time.

RNA-Seq has revolutionized the field of transcriptomics, enabling researchers to explore gene expression profiles, identify novel transcripts, discover
alternative splicing events, and investigate post-transcriptional modifications. It has applications in various areas of biological research, including
developmental biology, disease studies, drug discovery, and personalized medicine.


## Steps of RNA Seq before Sequencing

In brief, the basic steps of an RNA Seq experiment performed in the Lab are listed below:

1. RNA Extraction: Total RNA is extracted from the biological sample of interest, which could be cells, tissues, or even environmental samples.
2. RNA Fragmentation: The extracted RNA is fragmented into smaller pieces. This step helps to overcome the limitations of sequencing technologies, which have difficulty sequencing long RNA molecules.
3. cDNA Synthesis: The fragmented RNA is reverse-transcribed into complementary DNA (cDNA) using reverse transcriptase enzymes. These enzymes synthesize cDNA by using the RNA as a template.
4. Library Preparation: The cDNA is then prepared into a sequencing library. This involves adding specific adapters to the cDNA fragments, which are necessary for binding the fragments to the sequencing platform.
5. Sequencing: The prepared library is loaded onto a sequencing platform, such as Illumina or Ion Torrent sequencers. The fragments are then sequenced using high-throughput sequencing methods, generating millions of short sequence reads.

The process begins with the extraction of RNA from the sample of interest, which can be derived from various sources such as cells, tissues, or even single cells.
The extracted RNA is then converted into complementary DNA (cDNA) using reverse transcription. This step is necessary because most high-throughput sequencing
methods are more compatible with DNA rather than RNA.

Once the cDNA is obtained, it undergoes library preparation, where adapters are added to the ends of the cDNA fragments. These adapters contain specific sequences
that allow for the attachment of the cDNA fragments to a solid support, such as a flow cell or a glass slide.

Next, the cDNA fragments are amplified and undergo high-throughput sequencing. Several techniques can be used for sequencing, including Illumina sequencing, Ion
Torrent sequencing, and PacBio sequencing, among others. These sequencing platforms generate millions to billions of short sequences called reads.

If you want to read more about the basic conscepts of RNA Seq, you can always give a look to this [book](https://drive.google.com/file/d/10uR_guBPf75oF9lfAnu101llnnPIH2dx/view?usp=drive_link). Moreover, you can also ask google or chatGPT. 


Lets now dive to the dark side of RNA Seq.



## RNA Seq Data analysis

{% include image.html file="jekyll.png" url="http://jekyllrb.com" alt="Jekyll" caption="This is a sample caption" %}


After sequencing, the reads are processed to remove any artifacts or low-quality sequences. Then, they are aligned or mapped to a reference genome or
transcriptome to determine their origin and location within the genome. Alternatively, de novo assembly can be performed to reconstruct the transcriptome without
a reference genome.

The final step involves analyzing the mapped or assembled reads to quantify gene expression levels. This can be done by counting the number of reads that align to
specific genes or transcripts. Differential gene expression analysis can be performed to compare expression levels between different conditions or samples,
providing insights into the genes and biological pathways that are differentially regulated.


Read Mapping: The sequenced reads are aligned or mapped to a reference genome or transcriptome, depending on the experimental design. This step helps to determine
the origin and location of the RNA molecules.

Transcript Quantification: The number of reads mapped to each transcript is counted, providing information about the relative abundance of each RNA molecule in
the sample. This quantification can be used to compare gene expression levels between different samples or conditions.

Data Analysis: Bioinformatics analyses are performed on the generated data to interpret the results. This may involve identifying differentially expressed genes,
discovering novel transcripts, studying alternative splicing events, and performing functional enrichment analysis.


RNA-seq data analysis involves a series of computational steps to extract meaningful information from the generated sequencing data. Here is a general overview of
the RNA-seq data analysis pipeline:

1. Quality control (QC): The first step is to perform quality control on the raw sequencing reads. This involves assessing the quality scores of the bases,
  trimming low-quality bases or adapter sequences, and removing reads that are too short or contain sequencing artifacts.

2. Read alignment: The high-quality reads are then aligned or mapped to a reference genome or transcriptome using alignment algorithms such as HISAT2 or STAR.
  This step assigns the reads to specific genomic regions or transcripts.

3. Quantification: After alignment, the mapped reads are quantified to estimate the expression levels of genes or transcripts. This can be done using tools such
  as HTSeq, featureCounts. The output is typically a table that lists the number of reads assigned to each gene or transcript.

4. Normalization: To compare gene expression levels across samples, it is important to account for differences in sequencing depth and other technical factors. 
  Normalization methods such as DESeq2, edgeR, or TMM (trimmed mean of M-values) are commonly used to adjust for these factors and obtain more accurate expression 
  estimates.

5. Differential expression analysis: The normalized expression data can be used to identify genes that are differentially expressed between different experimental
  conditions or sample groups. Statistical methods like DESeq2, edgeR, or limma are employed to perform differential expression analysis, which helps identify
  genesthat show significant changes in expression levels.

6. Downstream Analysis:



Functional analysis: Once differentially expressed genes are identified, functional analysis can be performed to understand their biological significance. This
involves gene ontology (GO) analysis, pathway enrichment analysis, and other annotation tools to determine the biological processes, molecular functions, and
pathways associated with the differentially expressed genes.

Visualization: Visualization is an essential part of RNA-seq data analysis. Heatmaps, volcano plots, scatter plots, and other graphical representations help
visualize gene expression patterns, differential expression results, and clustering of samples or genes.

It is important to note that the specific tools and methods used in RNA-seq data analysis may vary depending on the research question, experimental design, and
available computational resources. Many software packages and workflows, such as the popular Galaxy platform, provide a user-friendly interface and preconfigured
pipelines to facilitate RNA-seq data analysis for researchers with limited bioinformatics expertise.



{% include links.html %}
