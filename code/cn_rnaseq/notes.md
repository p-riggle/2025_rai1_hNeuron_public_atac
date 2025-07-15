# 2023-03-05
Clustering of all DEGs -> ask what biological processes are enriched with each cluster. 
1. using tpm table, calcurate z score
2. then select only genes that showed differential expression between genotypes in any time point.  

# 2023-03-06
Heat map takes so much time and also PDF somehow does not save. I probably should reduce the number of genes significantly, figure out the saving then going back to the larger number of genes. 

# 2023-03-13
Clustering continues.

# 2023-03-23
Start doing pathway analysis of each cluster. 

# 2023-03-25
Compared pathways for each cluster.

# 2023-04-01

First column is gene names. 

Now, I would like to make a bar graph plotting the tpm values of a gene for each genotype (wt, het, or ko) grouped with days (d0 - d60). Each bar should denote the mean of all replicates. Along with each bar, I would like to put individual replicates' tpm values as small circles. And error bar should denote +/- standard error of mean from all replicates. Can you suggest an R script for such a plot for RAI1 gene? 

# 2023-04-25
I am going to make line plots of marker genes with Het and KO separetly.

# 2023-05-29
- Increasing threshold revealed opposite direction of gene misregulation between het and KO.
- Is it possible that larger fold changes mean direct target and smaller fold changes mean indirect target? -> I can plot frequency of fold changes. 

# 2023-05-31 SMSRF meeting
- In mouse data, what are the genes that are similar to human data?
- What are the evolutionary conserved RAI1-sensitive genes that are bound by RAI1 in mouse cortical neurons? -> the way to narrow down the drug target.  
- what are the genes in each cluster?
- what are the cellular compartment in each cluster? -> for Sexton lab collaboration
- what are the difference between weakly disregulated genes vs. strongly disregulated genes? 
- Is the genes similar between time points?
- Can we detect developmental acceleration in RAI1-KO mice? 

# 2023-06-13 Datasets dscription
List of Rai1 RNA-seq dataset
1. RAI1 CRISPRa in hippocampal neuron culture (DIV7-10): 
    Data: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP401688&o=acc_s%3Aa
    Paper: https://pubmed.ncbi.nlm.nih.gov/36410433/ (JBC 2023 Chang et al.)
2. RAI1 KO with Nestin Cre 3wk cortex: *they have striatum (two time points) and thalamus data with different cre. 
    Data: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP074487&o=acc_s%3Aa
    Paper: https://pubmed.ncbi.nlm.nih.gov/27693255/ (Neuron 2016 Huang et al.)
3. RAI1 KO mESCs differentating NPC (day 8) and neurons (day 12). 
    Data: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=ERP128947&o=acc_s%3Aa
    Paper: https://www.biorxiv.org/content/10.1101/2021.10.12.462922v1.full (likely, BioRxiv, 2023, Kaster et al.)
4. RAI1 RNAi in cortical neurons (DIV17). 2020 Cell Reports, Garay et al.
   
List of Rai1 ChIPseq dataset 
1. RAI1 ChIPseq cortical neurons (DIV17). 2020 Cell Reports, Garay et al. 
2. Flag-RAI1 ChIP-seq 8wk cortex:
    Data: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP074484&o=acc_s%3Aa
    Paper: https://pubmed.ncbi.nlm.nih.gov/27693255/ (Neuron 216 Huang et al.)

Strategy
1. Developmental accerelation. 
   A. mESC data 
   B. mouse Cortex data combined with encode data. 
2. Evolutionary conserved RAI1 targets 
   1. Select genes were disregulated in all the following data points.    
      1.  mESC-derived neurons at DIV12 (Kaster et al. 2023)
      2.  mouse cortex at 3wk (Huang et al. 2016)
      3.  hESC-derived neuron (this study)
   2. Of which pick common direct targets between in vitro and in vivo ChIP-seq data. 

But first do cellular component analysis of hESC mRNAseq data for Sexton collaboration. 

# 2023-06-15
- I noted DEG# discrepancy between DEG# plot and DEG# selected for clustering. DEG# selection for clustering was not correct because it used "|" between the samples.
- So fixed it. 
- no gap between k-means clusters in heatmap. Script fixed. But need to run again (it did not run as the server was unstable). Somehow k-means heatmap script freezes. Need to figure out why.

# 2023-06-16
- ran blueberry for Huang et al. 2016 data and Kaster et al. 2023 data.
- Strategy for picking genes that are disregulated in all the following data points.    
      1.  mESC-derived neurons at DIV12 (Kaster et al. 2023)
      2.  mouse cortex at 3wk (Huang et al. 2016)
      3.  hESC-derived neuron (this study)
1. DEG analysis for each dataset.
2. Aggregate DEGs from each dataset. Use the same threshold for all the datasets.
3. Look for overlapped genes. 
4. Change threshold for each dataset and repeat 1-3.

# 2023-06-18
- DEseq analysis of Huang and Kaster data sets

# 2025.06.27
- resumed the analysis!
- Organized UCSC genome browser sessions

# 2025.06.28
- Asked Satabdi to make line plot of Ciceri CRs 

# 2025.07.14
- started scRNAseq data analysis