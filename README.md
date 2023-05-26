# Pregnancy loss study

In this study we applied genome haplarithmisis to 94 miscarried products of conception (POCs) with normal parental and POC karyotypes. Utilizing parental DNA as well as POC extraembryonic mesoderm (EM) and chorionic villi (CV) DNA, representing embryonic and trophoblastic tissues, enabled characterization of the genomic landscape of both lineages. 

## 1.	Data conversion
SNP genotyping for the pregnancy loss study was performed in three runs. Raw genotyping data was converted into data that is compatible with haplarithmisis. 
Each run has its individual scripts to convert the raw genotyping data, however the steps and goal of these scripts generally are the same; 1. Extracting individual families based on unique PL code per family. 2. Converting column names to be haplarithmisis compatible. 

## 2.	Haplarithmisis
Haplarithmisis scripts are adapted from scripts developed by Masoud Zamani Esteki (PMID:25983246) and became publically available on publication by Masset et al. 2022 (PMID:35212381) https://github.com/GenomicsCoreLeuven/publications_JV/tree/main/GBS on Sep 13th 2021. 
The purpose here is to apply haplarithmisis to individual pregnancy loss families and determine the level of mosaicism, chimarism (maternal contamination). Haplarithmisis is a method that determines haplotypes, as well as the copy number and segregation origin of these haplotypes across the genome. The method advances and facilitates the detection of genuine DNA copy-number (or copy-neutral) aberrations and reveals their parental and mechanistic origin.

The produced haplarithm plots contain BAF, LogR, and haplotyping information on the fetal DNA samples (extraembryonic mesoderm and chorionic villi) and BAF and LogR of both parents. These plots anable the detection of copy number and copy neutral aberrations (CNV of >100kb, chromosomal aberrations, genome wide aberrations), parental (paternal or maternal) origin of the aberration, mechanistic (meiotic I, meitoic II, or mitotic) origin of the aberration, degree of mosaicism (range of 10%), and maternal cell contamination.

## 3.	CNV detection
The CNV detection script extracts genomic coordinates and aberrant segments based on LogR (copy-number) segments produced by haplarithmisis. The segments of logR are extracted by grouping consecutive genomic coordinates together that have identical logR values. The resulting table provides the genomic coordinates and size of any segmental aberration. The logR segmentation algorithm of haplarithmisis can be adjusted by changing the gamma values for LogR, in this study we applied the validated gamma of 14 by digtal droplet PCR that enabled us to detect large CNVs (>100 kb, PMID:31686035).  

## 4.	Mosaicism detection
### 4.1 Paternal/maternal haplotype segmentation 
The automated segmentation script extracts genomic coordinates and parental haplotyping data values leveraged and utilized by haplarithmisis. 
The extracted data consists of grouped P1, P2, M1, M2 haplarithm values, genomic coordinates, length and cumulative length. 
For mosaicism detection the genomic coordinates can be used for the 4.2 mean BAF calculations to indicate the region of interest or the P1, P2, M1, M2 values can be used directly for 4.3 mosaicism detection in the presence of maternal contamination.  

### 4.2 Mean BAF calculation
This script calculates the mean BAF values for the genomic coordinates of interest (obtained from 4.1 automated segmentation). It is required to indicate whether the aberration or region of interest shows 4 or 6 bands in the BAF plot. 
Resulting mean BAF values can be used in 4.3 to detect mosaicism degree in percentage. 

### 4.3 Mosacism detection 
These scripts compare either 1) mean BAF values (4 or 6 bands BAF), or 2) P1, P2, M1, M2 haplotyping values, to mosaicism tables from Conlin, et al., 2010 (PMID:20053666). 
In general, BAF values (1) are utilized for the mosaicism calculation. Only when maternal contamination is present, as identified in haplarithm plots, are the haplarithm values (2) utilized, to compensate for the disturbance in BAF values caused by maternal contamination. 

## 5.	Data analysis and visualization
Scripts for all main and extended data figures of the PL study. 
Main figure 1 A, B, C, D, main figure 2 B, and extended data figure 2 A, B, C, D.

## 6. Statistical analysis
Scripts containing the statistical analysis performed for this study.
