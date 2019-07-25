# `BCB420.2019.BioGRID`

#### (BioGRID data annotatation of human genes)

&nbsp;

## Installation

###### [Alison Wu] &lt;alison.wu@mail.utoronto.ca&gt;

## 1 About this package:


This package describe the workflow to download interaction datasets from [BioGRID](https://thebiogrid.org/),how to map the Entrez Gene IDs to HGNC symbols, and how to annotate the example gene set.


&nbsp;

#### In this project ...

```text
 --BCB420.2019.BioGrid/
   |__.gitignore
   |__.Rbuildignore
   |__BCB420.2019.BioGrid.Rproj
   |__DESCRIPTION
   |__dev/
      |__rptTwee.R
      |__toBrowser.R              
   |__inst/
      |__extdata/
         |entrez2sym.RData         # Entrez ID to HGNC symbol mapping tool
         |__xSetInt.tsv          # annotated example edges
      |__img/
         |__[...]                  # image sources for .md document
      |__scripts/
         |__recoverIDs.R           # utility to use biomaRt for ID mapping
   |__LICENSE
   |__NAMESPACE
   |__R/
      |__zzz.R
   |__README.md                    # this file

```

&nbsp;

----

## 2 BioGRID Data

BioGRID is a curated biological database of interactions. BioGRID interactions are recorded as the relationships between two genes or protiens which are supported by evidence code and pulication references.  These interactions include direct physical binding of two proteins, co-existence in a stable complex, and genetic interaction.

#### 2.1 Data semantics

Interactions that are recorded under the experimental systems (https://wiki.thebiogrid.org/doku.php/experimental_systems):

Genetic Interactions

1. Dosage Growth Defect – An overexpressed or increased dosage of one gene affects another gene that is mutated or deleted and causes a growth defect.
2. Dosage Lethality - An overexpressed or increased dosage of one gene affects another gene that is mutated or deleted and causes lethality in the stain. 
3.	Dosage Rescue – An overexpressed or increased dosage of one gene affects another gene that is mutated or deleted and reduces lethality in the stain.
4. Negative Genetic - A mutation or deletions in separate genes and each of genes causes a minimal phenotype and the same cell results in a severe fitness defect or lethality.   
5. Phenotypic Enhancement – A mutation/deletion/overexpression of one genes results in enhancement of any phenotype (other than lethality/growth defect) associated with mutation/deletion/overexpression of another gene, for example response to DNA damage or transcriptional output. 
6. Phenotypic Suppression – A mutation or deletion or overexpression of one gene results in suppression of any phenotype and results in mutation or deletion or overexpression of another gene.
7. Positive Genetic – A mutation or deletion in separate genes and one of the gene causes a minimal phenotype.  The same cell results in a less severe fitness defect than expected Synthetic Growth Defect - A mutation or deletion in separate genes and each of genes causes a minimal phenotype and the same cell results in a significant growth defect.
8. Synthetic Haploinsufficiency – A mutation or deletion in separate genes and one of each gene is hemizygous and cause a minimal phenotype and cause lethality
9. Synthetic Lethality - A mutation or deletion in separate genes and each of genes causes a minimal phenotype. 
10. Synthetic Rescue - A mutation or decreased dosage of one gene affects another gene that is mutated or deleted and reduces lethality or growth defect in the stain. 

Physical Interactions

1. Affinity Capture-Luminescence - An interaction is inferred when a bait protein that has been tagged with luciferase, is enzymatically detected in immunoprecipitates of the prey protein as light emission. 
2. Affinity Capture-MS - An interaction is inferred when the associated interaction partner is identified by mass spectrometric methods and a bait protein is affinity captured from cell extracts by polyclonal antibody or epitope tag. 
3. Affinity Capture-RNA - An interaction is inferred when a bait protein is captured from cell extracts by polyclonal antibody or epitope tag.
4. Affinity Capture-Western - An interaction is inferred when a bait protein is affinity captured from cell extracts by either polyclonal antibody or epitope tag and the associated interaction partner is identified by Western blot with a second epitope tag or specific polyclonal antibody.
5. Biochemical Activity - An interaction is inferred from the biochemical effect of two proteins. 
6. Co-crystal Structure - An interaction is inferred from the atomic level by X-ray crystallography. 
7. Co-fractionation - An interaction is inferred from the presence of multiple protein subunits in a partially purified protein preparation. 
8. Co-localization - An interaction is inferred from two proteins that co-localize in the cell by indirect immunofluorescence only when in addition, if one gene is deleted, the other protein becomes mislocalized. 
9. Co-purification - An interaction is inferred from the identification of two or more protein subunits in a purified protein complex. 
10. Far Western - An interaction is inferred when a bait protein is immobilized on a membrane and a prey protein that is incubated with the membrane localizes to the same membrane position as the bait protein. 
11. FRET - An interaction is inferred when close proximity of interaction partners is detected by fluorescence resonance energy transfer.
12. PCA - An interaction is inferred which a bait protein is expressed as a fusion to either an N- or C- terminal peptide fragment of a reporter protein and a prey protein is expressed as a fusion to the complementary C- or N- terminal fragment.
13. Protein-peptide - An interaction is inferred between a protein and a peptide derived from an interaction partner. 
14. Protein-RNA - An interaction is inferred between a protein and an RNA. 
15. Proximity Label-MS - An interaction is inferred when a bait-enzyme fusion protein selectively modifies a vicinal protein with a diffusible reactive product, followed by affinity capture of the modified protein and identification by mass spectrometric methods.
16. Reconstituted Complex - An interaction is inferred between proteins which include proteins in recombinant form or proteins isolated directly from cells with recombinant or purified bait. 
17. Two-hybrid - An interaction is inferred when a bait protein is expressed as a DNA binding domain fusion and when a prey protein is expressed as a transcriptional activation domain  fusion.




&nbsp;

## 3 Data download and cleanup

To download the source data from BioGRID :

1. Navigate to the [**BioGRID** database](https://thebiogrid.org/) 
2. Click "downloads" in the left corner, click "Current Release" and choose "BIOGRID-3.5.169" or follow the link [here](https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-3.5.174/)
3. Download the following data file: BIOGRID-ORGANISM-3.5.174.tab2.zip (71.19 MB)
4. Unzip the file
5. Choose the file `BIOGRID-ORGANISM-Homo_sapiens-3.5.174.tab2.txt` (130 MB) and place it into the your local working directory which is called `data`.

&nbsp;

## 4 Mapping Entrez IDs to HGNC symbols

BioGrid genes are Entrez IDs. These IDs can be mapped to HGNC symbols.

#### Preparations: packages, functions, files
Install a few required packages to begin:

1. **`readr`** allows the program to read large datasets.
```R
if (! requireNamespace("readr")) {
  install.packages("readr")
}
```

2. **`biomaRt`** is a Bioconductor package that implements the RESTful API of biomart,
the annotation framwork for model organism genomes at the EBI. It is a Bioconductor package, and as such it needs to be loaded via the **`BiocManager`**,
&nbsp;

```R
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
```

3. **`igraph`** is to compute statistics graphs and plots.
&nbsp;

```R
if (! requireNamespace("igraph")) {
  install.packages("igraph")
}
```

Fetch the HGNC reference data from GitHub(Steipe, 2019).

```R
myURL <- paste0("https://github.com/hyginn/",
                "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
load(url(myURL)) 
```


&nbsp;

#### 4.1 Step one: which IDs do we have to map?
&nbsp;

```R
  # Read the interaction data
  tmp <-readr::read_delim(file.path("../data","BIOGRID-ORGANISM-Homo_sapiens-3.5.174.tab2.txt"),
                          delim = "\t",
                          skip = 1, col_names = c("BioGRID Interaction ID", "Entrez_ID1", "Entrez_ID2"))  #468,058 rows
                          
  
  head(tmp)
  #    A tibble: 6 x 3
  #     `BioGRID Interaction ID` Entrez_ID1 Entrez_ID2
  #     1             <dbl>      <dbl>      <dbl>
  #     2              103       6416       2318
  #     3              117      84665         88
  #     4              183         90       2339
  #     5              278       2624       5371
  #     6              418       6118       6774
  #     7              586        375      23163                        
                          

  # The number of unique IDs we have to map:
  uENTREZ <- unique(c(tmp$Entrez_ID1, tmp$Entrez_ID2)) #24029 IDs
  
  ```

&nbsp;

#### 4.2  Step two: mapping via biomaRt

We first map Entrez IDs to HGNC symbols in bulk using biomaRt. And then we can look for the remaining IDs that we could not map via UniProt IDs (from HGNC reference data).


&nbsp;

###### 4.2.1  Constructing an ID-mapping tool

```R
# Map ENSP to HGNC symbols: open a "Mart" object ..
  myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

  tmp <- biomaRt::getBM(filters = "entrezgene_id",
                       attributes = c("entrezgene_id",
                                      "hgnc_symbol"),
                       values = uENTREZ,
                       mart = myMart)

  head(tmp)
  #      entrezgene_id hgnc_symbol
  # 1    10036      CHAF1A
  # 2    10048      RANBP9
  # 3    10059       DNM1L
  # 4    10084       PQBP1
  # 5    10128      LRPPRC
  # 6    10138        YAF2

  nrow(tmp)  #  16853 HGNC symbols have been retrieved for the 24029 Entrez IDs
  
```
&nbsp;

There are potential problems with the data that biomart returns

There might be duplicates (more than one value returned).
```R
sum(duplicated(tmp$entrezgene_id)) #87 duplicates
```

Or nothing returned for one Entrez ID
```R

  sum(! (uENTREZ) %in% tmp$entrezgene_id)  # 7263
```

&nbsp;

To fix the duplicates:

```R
  dupEntrez <- tmp$entrezgene_id[duplicated(tmp$entrezgene_id)]
  tmp[tmp$entrezgene_id %in% dupEntrez, ]
          entrezgene_id     hgnc_symbol
  # 366         6606            SMN2
  # 367         6606            SMN1
  # 465         9083            BPY2
  # 466         9083           BPY2B
  # 627        25788           BPY2C
  # 628        25788            FSBP
  # 935         8878          RAD54B
  # 936         8878          SQSTM1
  # 2046       23370           
  # 2047       23370        ARHGEF18
  # 2648        5414          
  # 2649        5414           SEPT4
  # 2657       54921           CHTF8
  # 2658       54921           DERPC
  # 2920        1124            CHN2
  
  

  #Remove the symbols with space
  tmp[tmp$hgnc_symbol %in% c(""), ]
  tmp <- tmp[ ! (tmp$hgnc_symbol %in% c("")), ]
  
  #Check the UniprotID for the duplicates one
  #SMN2 and SMN1 has the same UniprotID Q16637 - so assign arbitarily to either one
  tmp <- tmp[ ! (tmp$hgnc_symbol %in% c("SMN2","FSBP", "HLA-DQA2", "LINC02210-CRHR1", "DEFB4B", "DEFB103A", "CCL3L1", "CRYAA2", "IQCJ-SCHIP1", "TBC1D3L", "SERF1B", "DDT", "STAG3L3", "ZNF468", "KIR2DS5", "FNTB", "CHURC1-FNTB", "GTF2H2", "LGALS7B", "PLEKHG7", "FAM187A", "CBWD3", "BOLA2B", "RPP21", "TRIM39-RPP21", "LGALS7", "GUSBP1", "LINC00680","ARL17B", "USP17L2", "GATD3A", "RNA5S3", "PRR4", "GAGE12F", "GAGE12G")), ]
  
  #Check duplicates
  any(duplicated(tmp$entrezgene_id)) #FALSE
  ```
  
  &nbsp;
  
  Define the mapping tool after this cleanup
  
  &nbsp;
  
  ```R
  entrez2sym <- tmp$hgnc_symbol
  names(entrez2sym) <- tmp$entrezgene_id

  head(entrez2sym)
  # 10036    10048    10059    10084    10128    10138
  #"CHAF1A" "RANBP9"  "DNM1L"  "PQBP1" "LRPPRC"   "YAF2"
  ```
  &nbsp;
  
  ###### 4.2.2  Cleanup and validation of `entrez2sym`
  
  ```R
  sel <- ! (uENTREZ %in% names(entrez2sym))
 x <- rep(NA, sum( sel))
 names(x) <- uENTREZ[ sel ]

 any(duplicated(c(names(x), names(entrez2sym))))  #FALSE
entrez2sym <- c(entrez2sym, x)
 all(uENTREZ %in% names(entrez2sym)) #TRUE

-check if there are empty strings
sel <- which(entrez2sym == "")
 length(sel) #0

```
&nbsp;

Add out dated symbols

```R

  sel <- ( ! (entrez2sym %in% HGNC$sym)) & ( ! (is.na(entrez2sym)))
  length(        entrez2sym[ sel ] ) #101
  length( unique(entrez2sym[ sel ])) #101
  unkSym <- data.frame(unk = entrez2sym[ sel ],
                       new = NA,
                       stringsAsFactors = FALSE)

  View(unkSym)
  
  #add the outdated symbols(which are usually stored in the prev column in HGNC to entrez2sym)
  for (i in seq_len(nrow(unkSym))) {
      iPrev <- grep(unkSym$unk[i], HGNC$prev)[1]
      if (length(iPrev) == 1) {
          unkSym$new[i] <- HGNC$sym[iPrev]
          }
 }
 ```
 
 4.3 Final Validation
 
 Final validation for our mapping tool
 
 
```R
  all(uENTREZ %in% names(entrez2sym)) #TRUE, all the Entrez IDs are mapped
  
  sum(! is.na(entrez2sym)) #16671 symbols are found
  
  sum(! is.na(entrez2sym)) * 100 / length(entrez2sym) #69.34%
  
  #Save map
  save(entrez2sym, file = file.path("inst", "extdata", "entrez2sym.RData"))
  
  #Load file from RStudio project
  load(file = file.path("inst", "extdata", "entrez2sym.RData"))
  ```
  

&nbsp;

# 5 Data Statistics

Find the genes that are mostly show up in the both physical and genetic interactions (entrez2sym data).

```R
  HGNCsym = setNames(entrez2sym, NULL)
  num <- table(HGNCsym)
  num <- sort(num, decreasing=TRUE)
  barplot(num[1:10], main="Number of HGNC Symbol in entrez2sym")
```

![](./inst/img/num_of_HGNCsys.svg?sanitize=true "Number of HGNC Symbol in entrez2sym")


&nbsp;

# 6 Annotating gene sets with Example Data Set

We can now annotate the gene sets with BioGrid dta using our mapping tool. 



&nbsp;

```R
  tmp <-readr::read_delim(file.path("../data", "BIOGRID-ORGANISM-Homo_sapiens-3.5.174.tab2.txt"),
                          delim = "\t",
skip = 1, col_names = c("BioGRID Interaction ID", "Entrez_ID1", "Entrez_ID2", "ID1", "ID2", "sys_name1",   "sys_name2", "sym1", "sym2", "alias1", "alias2", "exp_name", "interaction_type", "author", "pubID", "org_ID1", "org_ID2", "interaction_throughput", "quant_score", "mod", "pheno", "qual", "tags", "ext"))

  #the number of interaction that is genetic
  nrow(tmp[tmp$interaction_type == "genetic",]) #7807
 
  #the number of interaction that is physical
  nrow(tmp[tmp$interaction_type == "physical",]) #489662
 
  #Finally we map Entrez IDs to HGNC sysmbols
  tmp$Entrez_ID1 <- entrez2sym[tmp$Entrez_ID1]
  tmp$Entrez_ID2 <- entrez2sym[tmp$Entrez_ID2]

  VALIDATE:
  sum(is.na(tmp$Entrez_ID1)) #178039
  sum(is.na(tmp$Entrez_ID2)) #180131
  # we remove edges in which either one or the other node is NA to
  # create our final data:
  BioGrid <- tmp[( ! is.na(tmp$Entrez_ID1)) & ( ! is.na(tmp$Entrez_ID2)), ]
  BioGridInt <- tmp[( ! is.na(tmp$Entrez_ID1)) & ( ! is.na(tmp$Entrez_ID2)), ]
  #Done
  save(BioGridInt, file = file.path(".", "data", "BioGridInt.RData"))
 
 
 #Network Graph
 t <- data.frame(Num_GeneticInteractions = num_genetic, Percentage=round(num_genetic / total * 100, 2), Num_PhysicalInteractions= num_phy,Percentage=round(num_phy / total * 100, 2),stringsAsFactors = FALSE)
  head(t)
  
  ```
  ![](./inst/img/interaction_percentage.svg?sanitize=true "Interaction Count and Percentage")

 ```R
 #the distribution
 deg <-  table(c(BioGridInt$Entrez_ID1, BioGridInt$Entrez_ID2))
 
  hist(deg, breaks=50,
      xlim = c(0, 1400),
      col = "#3fafb388",
      main = "Degree distribution",
      xlab = "degree (undirected graph)",
      ylab = "Counts")
 ```
      
  ![](./inst/img/degree_dist.svg?sanitize=true "Degree Distribution")

&nbsp;



To conclude, we annotate the example gene set, validate the annotation, and store the data.

&nbsp;

```R

  # The specification of the sample set is copy-paste from the 
  # BCB420 resources project.

  xSet <- c("AMBRA1", "ATG14", "ATP2A1", "ATP2A2", "ATP2A3", "BECN1", "BECN2",
          "BIRC6", "BLOC1S1", "BLOC1S2", "BORCS5", "BORCS6", "BORCS7",
          "BORCS8", "CACNA1A", "CALCOCO2", "CTTN", "DCTN1", "EPG5", "GABARAP",
          "GABARAPL1", "GABARAPL2", "HDAC6", "HSPB8", "INPP5E", "IRGM",
          "KXD1", "LAMP1", "LAMP2", "LAMP3", "LAMP5", "MAP1LC3A", "MAP1LC3B",
          "MAP1LC3C", "MGRN1", "MYO1C", "MYO6", "NAPA", "NSF", "OPTN",
          "OSBPL1A", "PI4K2A", "PIK3C3", "PLEKHM1", "PSEN1", "RAB20", "RAB21",
          "RAB29", "RAB34", "RAB39A", "RAB7A", "RAB7B", "RPTOR", "RUBCN",
          "RUBCNL", "SNAP29", "SNAP47", "SNAPIN", "SPG11", "STX17", "STX6",
          "SYT7", "TARDBP", "TFEB", "TGM2", "TIFA", "TMEM175", "TOM1",
          "TPCN1", "TPCN2", "TPPP", "TXNIP", "UVRAG", "VAMP3", "VAMP7",
          "VAMP8", "VAPA", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39",
          "VPS41", "VTI1B", "YKT6")


  # Example genes are not among the known nodes?
  x <- which( ! (xSet %in% c(BioGridInt$Entrez_ID1, BioGridInt$Entrez_ID2)))
  cat(sprintf("\t%s\t(%s)\n", HGNC[xSet[x], "sym"], HGNC[xSet[x], "name"]))


# ATP2A1	(ATPase sarcoplasmic/endoplasmic reticulum Ca2+ transporting 1)
# 	ATP2A2	(ATPase sarcoplasmic/endoplasmic reticulum Ca2+ transporting 2)
# BECN2	(beclin 2)
# BIRC6	(baculoviral IAP repeat containing 6)
# BLOC1S2	(biogenesis of lysosomal organelles complex 1 subunit 2)
# BORCS5	(BLOC-1 related complex subunit 5)
# BORCS7	(BLOC-1 related complex subunit 7)
# BORCS8	(BLOC-1 related complex subunit 8)
# CTTN	(cortactin)
# EPG5	(ectopic P-granules autophagy protein 5 homolog)
# GABARAPL2	(GABA type A receptor associated protein like 2)
# INPP5E	(inositol polyphosphate-5-phosphatase E)
# IRGM	(immunity related GTPase M)
# LAMP5	(lysosomal associated membrane protein family member 5)
# MAP1LC3B	(microtubule associated protein 1 light chain 3 beta)
# MGRN1	(mahogunin ring finger 1)
# OPTN	(optineurin)
# OSBPL1A	(oxysterol binding protein like 1A)
# PI4K2A	(phosphatidylinositol 4-kinase type 2 alpha)
# PIK3C3	(phosphatidylinositol 3-kinase catalytic subunit type 3)
# RAB20	(RAB20, member RAS oncogene family)
# RAB29	(RAB29, member RAS oncogene family)
# RAB34	(RAB34, member RAS oncogene family)
# RAB7A	(RAB7A, member RAS oncogene family)
# RPTOR	(regulatory associated protein of MTOR complex 1)
# RUBCNL	(rubicon like autophagy enhancer)
# SNAP29	(synaptosome associated protein 29)
# SNAP47	(synaptosome associated protein 47)
# SPG11	(SPG11, spatacsin vesicle trafficking associated)
# STX17	(syntaxin 17)
# STX6	(syntaxin 6)
# TFEB	(transcription factor EB)
# TGM2	(transglutaminase 2)
# TOM1	(target of myb1 membrane trafficking protein)
# TPCN1	(two pore segment channel 1)
# TPCN2	(two pore segment channel 2)
# TPPP	(tubulin polymerization promoting protein)
# VPS11	(VPS11, CORVET/HOPS core subunit)
# YKT6	(YKT6 v-SNARE homolog)

  #select the gene  that are in the example set:
  sel <- (BioGridInt$Entrez_ID1 %in% xSet) & (BioGridInt$Entrez_ID2 %in% xSet)
  xSetInt <- BioGridInt[sel, c("Entrez_ID1", "Entrez_ID2")]

  #Statitics:
  nrow(xSetInt) #36

  # Save the annotated set
  writeLines(c("Entrez_ID1\tb",
                            sprintf("%s\t%s", xSetInt$Entrez_ID1, xSetInt$Entrez_ID2)),
                        con = "xSetInt.tsv")
  # The data set can be read back
  myXset <- read.delim(file.path("inst", "extdata", "xSetInt.tsv"),
                      stringsAsFactors = FALSE)


  myXset <- read.delim(system.file("extdata",
                                  "xSetInt.tsv",
                                  package = "BCB420.2019.BioGrid"),
                      stringsAsFactors = FALSE)
 ```

&nbsp;
