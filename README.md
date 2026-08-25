# IFN646 Project

## Overview

In the project you will work as a group to complete a multi-stage analysis of real-world biomedical data. The work will involve some concepts we directly cover in the unit, and others you will learn while doing the project.

You are working as part of a group, typically of four students. Students who do not form their own group are placed in groups by the teaching team.

## Timeline

- Week 4: group formation starts

- Week 5: project details are discussed, and groups are finalised. The tutorial is used to discuss the assessment.

- Week 7: prac slots are used for project support

- Week 10: tutorial and prac slots are used for project support

- Week 12: tutorial slot is used for support. Project reports are due

- Week 13: oral presentations take place (during the timetabled activities)

## Data

- RNA-seq raw counts

- Human Genome assembly GRCh38: [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/)

  + Genomic data (GCF_000001405.26_GRCh38_genomic.fna.gz, 902MB compressed)
  
  + Genome annotation (GCF_000001405.26_GRCh38_genomic.gff.gz, 24MB compressed)

- Variant information from the '1000 genomes' project

  + You can view the whole collection on the project [FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/)
  
  + We are providing a preprocessed subset for you on this page:

    - Variants from 7 individuals on chromosome 1: sampled_chr1.vcf

    - Variants from the same individuals on chromosome 2: sampled_chr2.vcf

    - Each file is just over 110MB.

Note that not all files are needed for all steps in the project.

## Tasks

### Context and structure

In this project you will work as a group to complete a multi-stage analysis of real-world biomedical data. You will collectively play the role of a data scientist collaborating on a biomedical research project.

The data provided above is completely real, and all the tasks you will complete are authentic, but their selection and arrangement have been designed to support your learning and our assessment of your learning.

The overall objective of the collaboration is to develop novel drugs against COVID-19. Of course, this is a 'pretend' collaboration, and we are not developing drugs as part of the project. The underlying hypothesis that frames the project is that:

1. COVID-19 patients suffer from a deregulation in gene expression

2. CRISPR can be used to regulate the expression of specific genes of interest

3. If we can restore the expression of specific genes deregulated in COVID-19, it will lead to a faster and more complete recovery

The first part of the hypothesis is true. You will confirm this for yourself in task 1.

The second part is also correct. It is possible to regulate transcription using CRISPRa (CRISPR activation) and CRISPRi (CRISPR interference), as discussed [here](https://genome.cshlp.org/content/31/11/2120). It has been used in other contexts, such as [cancer](https://www.nature.com/articles/s41590-019-0500-4) and [neurodegenerative diseases](https://www.nature.com/articles/srep28420).

What is 'made up' to provide a structure to the project is the third part.

### Task 1

Your starting point will be an RNA-seq experiment, and you will have to identify genes with a significantly difference in expression in COVID-19. You want to be certain, so you will use the consensus between at least two tools.

You will need to decide:

- Which tools to consider: Other than DESeq2, what other tools are available? What are their strengths and limitations? Which criteria should be used to select the tool?
    
- How to run them
    
- How to combine the results: Each tool will provide its set of results. How to combine them? What strategies can you use to achieve this?

You are free to use DESeq2 as one of the tools in this task, but you do not have to. It is entirely up to you.

### Task 2

Next, you will select one of the genes you have identified. For this gene, you will have to design suitable guide RNA (gRNA) sequences for CRISPR (either CRISPRi or CRISPRa, depending on the direction of the change of expression)

You will need to extract the DNA sequence for the gene you want to target. You can do this using the genome assembly data shared above and a similar approach to the week 1 workshop. You are also free to use any other suitable resource, such as NCBI GeneLinks to an external site.

Different CRISPRi/CRISPRa systems may have different constraints in terms of the exact location to target. To simplify task 2, your objective is to identify suitable gRNAs that target any region in the first exon of the gene you have selected. Please note that we are talking about the first exon relative to the start of the gene: you need to pay attention to directionality.

- On-target evaluation: How can we maximise the chances of obtaining the desired edits? In principle, editing should be easy: (1) Take your favorite gene, (2) Identify a suitable sequence, (3) Construct the guide RNA, (4) Voila

- Off-target evaluation: How can we minimise the risk of off-target modifications? 

    + Be careful with those scissors!

        + We want to target a given gene, but nothing else
        
        + We need a high degree of confidence on where we are cutting
        
        + It is essential for sequence of the guide RNA to be unique
        
        + It is not enough
        
        + We need to be aware of potential off-target modifications
    
    + Partial matches are dangerous
    
        + Off-target sequences with a small number of mismatches to the gRNA sequence can lead to problems.
        
        + Some tools report, for any candidate gRNA, the number of off-target sites with 1, or 2, or 3, or 4 mismatches, but not all off-target site are equally dangerous.
    
    + A more refined assessment
    
        + The number of mismatches is important, but their position matter too
        
        + For a given guide RNA, we can extract all potential off-target sites with at most 4 mismatches
        
        + For each pair of the guide and an off-target site, we can calculate a score based on the number of mismatches, their position, and the distance between them
        
        + We can then combine all these local scores into a global score for the guide
        
        + This approach is known as the MIT score.
        
        + The type of mismatch (e.g. G→C or G→T) also plays a role
        
        + As before, we can extract highly-similar off-target sites, calculate a local score for each pair of the guide and an off-target site, and combine these local scores into a global score for the guide
        
        + This approach is known as the CFD score. Other approaches exist too.

A number of gRNA design tools are available. You only need to use one tool, and you are free to select which one to use, but you must justify your choice. Your choice of tool may be guided by what you also have to complete in task 3.

For task 2, we ask to pay attention to the predicted efficiency of the guides (i.e., will they work?) as well as their specificity (i.e., will they lead to any off-target effects?).

### Task 3

Finally, you will work with real human genomes and refine the analysis completed in task 2.

We all have slightly different genomes, and it will impact the suitability of the gRNAs you have previously identified.

In particular, you will need to assess whether:

- some guides may no longer be present in some of these genomes
    
- some guides may have a different risk of off-target modifications in some of these genomes

The genomes are made available as a set of VCF files that capture the differences between these genomes and the reference genome. You need to decide how you use that information to address these two questions.

## Submission

There are two separate submission links for your project:

- You will submit as group a final report, worth 40 marks
    
- You will individually submit a review/reflection, worth 15 marks

## Structure and Templates

As different groups may be using different tools to produce the report (Word, LaTeX, notebooks, etc.), we are not providing a template, but we provide a recommended structure.

There is no page limit for the report, but remember that we mark the quality of the content, not the quantity. If something can be explained concisely, please do.

Similarly, there is no page limit for the reflection, but we would expect most submissions to be approximately 1-2 pages.
