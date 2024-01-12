---
title: Gene ID Conversion
author: ''
date: '2024-01-07'
slug: []
categories: []
tags: []
excerpt: ""
---

Various ways to convert ENTREZ ID to gene symbols and vice versa


```r
library(tidyverse)
library(DOSE)
```

# load dataset


```r
# load a list of ENTREZ ID's
data(geneList)
geneList %>% head
```

```
##     4312     8318    10874    55143    55388      991 
## 4.572613 4.514594 4.418218 4.144075 3.876258 3.677857
```

```r
length(names(geneList))
```

```
## [1] 12495
```

# biomart


```r
library(biomaRt)
```


```r
biomaRt::listEnsemblArchives()
```

```
##              name     date                                 url version
## 1  Ensembl GRCh37 Feb 2014          https://grch37.ensembl.org  GRCh37
## 2     Ensembl 110 Jul 2023 https://jul2023.archive.ensembl.org     110
## 3     Ensembl 109 Feb 2023 https://feb2023.archive.ensembl.org     109
## 4     Ensembl 108 Oct 2022 https://oct2022.archive.ensembl.org     108
## 5     Ensembl 107 Jul 2022 https://jul2022.archive.ensembl.org     107
## 6     Ensembl 106 Apr 2022 https://apr2022.archive.ensembl.org     106
## 7     Ensembl 105 Dec 2021 https://dec2021.archive.ensembl.org     105
## 8     Ensembl 104 May 2021 https://may2021.archive.ensembl.org     104
## 9     Ensembl 103 Feb 2021 https://feb2021.archive.ensembl.org     103
## 10    Ensembl 102 Nov 2020 https://nov2020.archive.ensembl.org     102
## 11    Ensembl 101 Aug 2020 https://aug2020.archive.ensembl.org     101
## 12    Ensembl 100 Apr 2020 https://apr2020.archive.ensembl.org     100
## 13     Ensembl 99 Jan 2020 https://jan2020.archive.ensembl.org      99
## 14     Ensembl 98 Sep 2019 https://sep2019.archive.ensembl.org      98
## 15     Ensembl 97 Jul 2019 https://jul2019.archive.ensembl.org      97
## 16     Ensembl 96 Apr 2019 https://apr2019.archive.ensembl.org      96
## 17     Ensembl 95 Jan 2019 https://jan2019.archive.ensembl.org      95
## 18     Ensembl 94 Oct 2018 https://oct2018.archive.ensembl.org      94
## 19     Ensembl 93 Jul 2018 https://jul2018.archive.ensembl.org      93
## 20     Ensembl 80 May 2015 https://may2015.archive.ensembl.org      80
## 21     Ensembl 77 Oct 2014 https://oct2014.archive.ensembl.org      77
## 22     Ensembl 75 Feb 2014 https://feb2014.archive.ensembl.org      75
## 23     Ensembl 54 May 2009 https://may2009.archive.ensembl.org      54
##    current_release
## 1                 
## 2                *
## 3                 
## 4                 
## 5                 
## 6                 
## 7                 
## 8                 
## 9                 
## 10                
## 11                
## 12                
## 13                
## 14                
## 15                
## 16                
## 17                
## 18                
## 19                
## 20                
## 21                
## 22                
## 23
```

```r
listMarts()
```

```
##                biomart                version
## 1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 110
## 2   ENSEMBL_MART_MOUSE      Mouse strains 110
## 3     ENSEMBL_MART_SNP  Ensembl Variation 110
## 4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 110
```

```r
mart_aug2020 <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                         dataset = "hsapiens_gene_ensembl",
                         host = "https://aug2020.archive.ensembl.org")

head(listAttributes(mart_aug2020), n = 10)
```

```
##                             name                  description         page
## 1                ensembl_gene_id               Gene stable ID feature_page
## 2        ensembl_gene_id_version       Gene stable ID version feature_page
## 3          ensembl_transcript_id         Transcript stable ID feature_page
## 4  ensembl_transcript_id_version Transcript stable ID version feature_page
## 5             ensembl_peptide_id            Protein stable ID feature_page
## 6     ensembl_peptide_id_version    Protein stable ID version feature_page
## 7                ensembl_exon_id               Exon stable ID feature_page
## 8                    description             Gene description feature_page
## 9                chromosome_name     Chromosome/scaffold name feature_page
## 10                start_position              Gene start (bp) feature_page
```

```r
# find the entrez gene id name
grep( "entrez", listAttributes(mart_aug2020)$name, value = T)
```

```
## [1] "entrezgene_trans_name"  "entrezgene_description" "entrezgene_accession"  
## [4] "entrezgene_id"
```


```r
t2g_aug2020 <- biomaRt::getBM(attributes = c("entrezgene_id","external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version", "gene_biotype", "transcript_biotype","chromosome_name", "band", "transcript_length", "start_position", "end_position","transcription_start_site", "strand", "refseq_mrna","refseq_ncrna"), mart = mart_aug2020)

# rename
t2g_aug2020 <- dplyr::rename(t2g_aug2020,  entrez_gene = entrezgene_id, ext_gene = external_gene_name, ens_gene = ensembl_gene_id, ens_gene_ver = ensembl_gene_id_version)
```


```r
# biomart genesymbols
index <- match(names(geneList), t2g_aug2020$entrez_gene)
biomart_gs <- t2g_aug2020[index,] %>% .$ext_gene
# wasn't able to use this because the NA's from match is ignored
# biomart_gs <- t2g_aug2020 %>% slice(match(names(geneList), entrez_gene)) %>% .$ext_gene
```

# annotationhub


```r
library(AnnotationHub)
```


```r
# these codes are from https://github.com/hbctraining/scRNA-seq/blob/master/lessons/mitoRatio.md

# Connect to AnnotationHub
# the snapshot is saved at /home/aahn/.cache/R/AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

#Next, we acquire the latest annotation files from this Ensembl database.
#We can first check which annotation versions are available:
# Check versions of databases available
ahDb %>% 
  mcols()

# Since we want the most recent, we will return the AnnotationHub ID for this database:
# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

#Finally, we can use the AnnotationHub connection to download the appropriate Ensembl database, which should be version GRCh38.92.

# Download the appropriate Ensembldb database
# this took ages so will need to save the rds 
edb <- ah[[id]]

# this did not work ... i guess you cant save a ENsDb object? 
#saveRDS(edb, "/researchers/antonio.ahn/resources/R_resources/annotationhub/data/edb_AH104864.rds")

#And to extract gene-level information we can use the Ensembldb function genes() to return a data frame of annotations.

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")       

# saved as these steps takes quite a while
saveRDS(annotations, "/researchers/antonio.ahn/resources/R_resources/annotationhub/data/annotations_AH104864.rds")
```


```r
annoHub <- readRDS("/researchers/antonio.ahn/resources/R_resources/annotationhub/data/annotations_AH104864.rds")
```


```r
index <- match(names(geneList), annoHub$entrezid)
annoHub_gs <- annoHub[index,] %>% .$symbol
```

# mapIds


```r
library(org.Hs.eg.db)
```


```r
mapIds_gs <- mapIds(org.Hs.eg.db, keys=names(geneList), keytype="ENTREZID", columns="SYMBOL",column="SYMBOL")
```

# summary

in this example, the mapID method gave the most non NA's


```r
summary_df <- tibble(entrezID = names(geneList), biomart_gs = biomart_gs, mapIds_gs, annoHub_gs)
```


```r
summary_df$biomart_gs %>% is.na %>% summary
```

```
##    Mode   FALSE    TRUE 
## logical   12263     232
```

```r
summary_df$mapIds_gs %>% is.na %>% summary
```

```
##    Mode   FALSE    TRUE 
## logical   12440      55
```

```r
summary_df$annoHub_gs %>% is.na %>% summary
```

```
##    Mode   FALSE    TRUE 
## logical   12214     281
```



