# Pregnancy_loss_study
## 1.	DataConversion
SNP genotyping for the pregnancy loss study was performed in three runs. Raw genotyping data was converted into data that is compatible with haplarithmisis. 
Each run has its individual scripts to convert the raw genotyping data, however the steps and goal of these scripts generally are the same; 1. Extracting individual families based on unique PL code per family. 2. Converting column names to be haplarithmisis compatible. 

## 2.	Haplarithmisis
Haplarithmisis scripts are adapted from scripts by Masoud Zamani Esteki. 
The purpose is to apply haplarithmisis to individual pregnancy loss families. Haplarithmisis is a method that determines haplotypes, as well as the copy number and segregation origin of these haplotypes across the genome. The method advances and facilitates the detection of genuine DNA copy-number (or copy-neutral) aberrations and reveals their parental and mechanistic origin.

The produced haplarithm plots contain BAF, LogR, and haplotyping information on the fetal DNA samples (Extraembryonic mesoderm and chorionic villi) and BAF and LogR of both parents. These plots anable the detection of copy number and copy neutral aberrations (CNV of >100kb, chromosomal aberrations, genome wide aberrations), parental (paternal or maternal) origin of the aberration, mechanistic (meiotic I, meitoic II, or mitotic) origin of the aberration, degree of mosaicism (range of 10%), and maternal cell contamination.

More on haplarithmisis: 
https://doi.org/10.1016/j.ajhg.2015.04.011
This paper is specifically on haplarithmisis on single cells however the pregnancy loss study utilizes bulk DNA samples. 

## 3.	CNV_Detection

## 4.	MosaicismDetection

## 5.	DataAnalysis_Plotting_Statistics

## Miscellaneous
