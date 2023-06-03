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

RNA Seq is a biology experiment used for quantification of gene expression. The goal of this method is to define which genes are transcribed and which are not in
a specific time and in a specific cell type (population). Quantification refers to not only determine which genes are transcribed, but also their abundance,
giving us the opportunity to study the biological processes affecting the abundance of different DNA fragments.

## Steps of RNA Seq before Sequencing

In brief, the basic steps (lets say wet lab steps) of an RNA Seq experiment are listed below:

1. RNA extraction: Extraction of all the cell's RNA.
2. RNA to cDNA (complementary DNA) transition via a process known as reverse transcription.
3. Sequence of cDNA.

After sequencing those cDNAs, we observe the abundances of DNA, and then we attempt to infer the original amounts of RNA in the cell. Despite the fact that this
task seems quite simple, in reality it is not. 

If you want to read more about the basic conscepts of RNA Seq, you can always give a look to this [book](https://drive.google.com/file/d/10uR_guBPf75oF9lfAnu101llnnPIH2dx/view?usp=drive_link). Moreover, you can also ask google or chatGPT. 


Lets now dive to the dark side of RNA Seq.



## RNA Seq Data analysis













## About Ruby Gems

Ruby has a number of plugins referred to as "gems." Just because you have Ruby doesn't mean you have all the necessary Ruby gems that your program needs to run. Gems provide additional functionality for Ruby programs. There are thousands of [Rubygems](https://rubygems.org/) available for you to use.

Bundler looks in a project's "Gemfile" (no file extension) to see which gems are required by the project. The Gemfile lists the source and then any gems, like this:

```
source "https://rubygems.org"

gem 'github-pages'
gem 'jekyll'
```

The source indicates the site where Bundler will retrieve the gems: [https://rubygems.org](https://rubygems.org).

The gems it retrieves are listed separately on each line.

Here no versions are specified. Sometimes gemfiles will specify the versions like this:

```
gem 'kramdown', '1.0'
```


See this [Stack Overflow post](http://stackoverflow.com/questions/5170547/what-does-tilde-greater-than-mean-in-ruby-gem-dependencies) for more details.


{% include links.html %}
