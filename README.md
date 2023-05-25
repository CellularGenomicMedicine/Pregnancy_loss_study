# Pregnancy loss study

In this study we applied genome haplarithmisis to 94 miscarried products of conception (POCs) with normal parental and POC karyotypes. Utilizing parental DNA as well as POC extraembryonic mesoderm (EM) and chorionic villi (CV) DNA, representing embryonic and trophoblastic tissues, enabled characterization of the genomic landscape of both lineages. In contrast to viable pregnancies where mosaic chromosomal abnormalities are often restricted to CV, such as confined placental mosaicism, we found that in PLs, the situation is reversed with a higher degree of mosaic chromosomal imbalances in EM rather than CV, indicating that the aberrations originate by the end of the first embryonic week, before blastocyst formation. Our results stress the critical importance of scrutinizing the full allelic architecture of genomic abnormalities in pregnancy loss to improve the clinical management and basic research of this devastating condition.

## 1.	Data conversion
SNP genotyping for the pregnancy loss study was performed in three runs. Raw genotyping data was converted into data that is compatible with haplarithmisis. 
Each run has its individual scripts to convert the raw genotyping data, however the steps and goal of these scripts generally are the same; 1. Extracting individual families based on unique PL code per family. 2. Converting column names to be haplarithmisis compatible. 

## 2.	Haplarithmisis
Haplarithmisis scripts are adapted from scripts by Masoud Zamani Esteki. 
The purpose is to apply haplarithmisis to individual pregnancy loss families. Haplarithmisis is a method that determines haplotypes, as well as the copy number and segregation origin of these haplotypes across the genome. The method advances and facilitates the detection of genuine DNA copy-number (or copy-neutral) aberrations and reveals their parental and mechanistic origin.

The produced haplarithm plots contain BAF, LogR, and haplotyping information on the fetal DNA samples (Extraembryonic mesoderm and chorionic villi) and BAF and LogR of both parents. These plots anable the detection of copy number and copy neutral aberrations (CNV of >100kb, chromosomal aberrations, genome wide aberrations), parental (paternal or maternal) origin of the aberration, mechanistic (meiotic I, meitoic II, or mitotic) origin of the aberration, degree of mosaicism (range of 10%), and maternal cell contamination.

More on haplarithmisis: 
https://doi.org/10.1016/j.ajhg.2015.04.011
This paper is specifically on haplarithmisis applied to single cells however the pregnancy loss study utilizes bulk DNA samples. 

## 3.	CNV detection
The CNV detection script extracts genomic coordinates and aberrant segments based on LogR (copy-number) segments produced by haplarithmisis. The segments of logR are extracted by grouping consecutive genomic coordinates together that have identical logR values. The resulting table provides the genomic coordinates and size of any segmental aberration. The logR segmentation algorithm of haplarithmisis can be adjusted by changes the gamma values for LogR, in the pregnancy loss study gamma of 14 is used. This value requires optimization for the specific samples that are being analyzed. 

## 4.	Mosaicism detection
### 4.1 Paternal/maternal haplotype segmentation 
The automated segmentation script extracts genomic coordinates and parental haplotyping data values leveraged and utilized by haplarithmisis. 
The extracted data consists of grouped P1, P2, M1, M2 haplotyping values, genomic coordinates, length and cumulative length. 
For mosaicism detection the genomic coordinates can be used for the 4.2 mean BAF calculations to indicate the region of interest or the P1, P2, M1, M2 values can be used directly for 4.3 mosaicism detection in the presence of maternal contamination.  

### 4.2 Mean BAF calculation
This script calculates the mean BAF values for the genomic coordinates of interest (obtained from 4.1 automated segmentation). It is required to indicate whether the aberration or region of interest shows 4 or 6 bands in the BAF plot. 
Resulting mean BAF values can be used in 4.3 to detect mosaicism degree in percentage. 

### 4.3 Mosacism detection 
These scripts compare either 1) mean BAF values (4 or 6 bands BAF), or 2) P1, P2, M1, M2 haplotyping values, to mosaicism tables from Conlin, et al., 2010. 
In general, BAF values (1) are utilized for the mosaicism calculation. Only when maternal contamination is present, as identified in haplarithm plots, are the haplotyping values (2) utilized, to compensate for the disturbance in BAF values caused by maternal contamination. 

## 5.	Data analysis, visualization and statistics

## Miscellaneous
