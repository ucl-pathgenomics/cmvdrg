
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cmvdrg

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/ucl-pathgenomics/cmvdrg.svg?branch=master)](https://travis-ci.org/ucl-pathgenomics/cmvdrg)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/)
![GitHub repo
size](https://img.shields.io/github/repo-size/ucl-pathgenomics/cmvdrg.svg)
![GitHub code size in
bytes](https://img.shields.io/github/languages/code-size/ucl-pathgenomics/cmvdrg.svg)
<!-- badges: end --> The goal of cmvdrg is to enable antiviral drug
Resistance Genotyping, for Human Cytomegalovirus sequencing data.

## Installation

You can install the current version from
[GitHub](https://github.com/ucl-pathgenomics/cmvdrg) with:

``` r
# install.packages("devtools")
devtools::install_github("ucl-pathgenomics/cmvdrg")
```

## Database

A manually curated csv file, containing the relationships between - Gene
- AA mutation - Drug susceptibility, values are EC50 fold changes to
susceptible strains. “Susceptible” or “Resistant” are present if data is
anecdotal. - The method for Resistance Phenotyping, including where
possible strian information used in marker transfer. i.e. AD169, TOWNE.
- Publication reference, DOI, web address & review paper reference if
used. - Where a multi-gene, multi AA mutation susceptibility profile is
available for a known drug, this is included too. - Administrative
information, last update to reference, updater, whether the row is used
in calculations or not. i.e. if a mutations is present in a review
paper, but is not traceable to given references this column is set to
“U” Unused, or “R” to Review.

## Web service

A user friendly Shiny Applications has been bundled with this package.
The same application is available over the internet here
<http://51.11.13.133:3838/cmvdrg/> where the terms of use are contained.

## Example

Load the package for usage

``` r
library("cmvdrg")
```

### Call resistance mutations from cmvdrg

``` r
# load some test data, bundled with the package
# a varscan tab file & fasta files are also bundled as examples.
my_sample = system.file("testdata", "A10.vcf", package = "cmvdrg")

# this calls any variants in the sample file, that are also present in the cmvdrg database.
data = call_resistance(infile = my_sample, all_mutations = F)
#> Loading required package: GenomicFeatures
#> Loading required package: BiocGenerics
#> Loading required package: parallel
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:parallel':
#> 
#>     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
#>     clusterExport, clusterMap, parApply, parCapply, parLapply,
#>     parLapplyLB, parRapply, parSapply, parSapplyLB
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
#>     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
#>     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
#>     union, unique, unsplit, which, which.max, which.min
#> Loading required package: S4Vectors
#> Loading required package: stats4
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:base':
#> 
#>     expand.grid
#> Loading required package: IRanges
#> Loading required package: GenomeInfoDb
#> Loading required package: GenomicRanges
#> Loading required package: AnnotationDbi
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 'select()' returned many:1 mapping between keys and columns

# view the resutant EC50 fold changes to Gangiclovir.
# if you have NGS data you have have informative mutational frequencies also.
print(data[,c("change", "freq", "Ganciclovir")])
#>        change   freq Ganciclovir
#> 1  UL54_D588N   5.2%           2
#> 2  UL54_D588N   5.2%         3.8
#> 3  UL54_N685S   100% Susceptible
#> 4  UL54_S655L   100% Susceptible
#> 5  UL97_C592G 10.77%         2.9
#> 6  UL97_C592G 10.77%           3
#> 7  UL97_C592G 10.77%           3
#> 8  UL97_H411Y  1.69%            
#> 9  UL97_H411Y  1.69%         0.5
#> 10 UL97_T409M 37.13%            
#> 11 UL97_T409M 37.13%         0.9

# view the first few lines, there are a lot of columns!
head(data)
#>       change    seqnames  start    end width strand         id   freq RefCount
#> 1 UL54_D588N NC_006273.2  80161  80161     1      - single run   5.2%      529
#> 2 UL54_D588N NC_006273.2  80161  80161     1      - single run   5.2%      529
#> 3 UL54_N685S NC_006273.2  79869  79869     1      - single run   100%        0
#> 4 UL54_S655L NC_006273.2  79959  79959     1      - single run   100%        0
#> 5 UL97_C592G NC_006273.2 143571 143571     1      + single run 10.77%      174
#> 6 UL97_C592G NC_006273.2 143571 143571     1      + single run 10.77%      174
#>   VarCount VarAllele varAllele CDSLOC.start CDSLOC.end CDSLOC.width PROTEINLOC
#> 1       29         T         A         1762       1762            1        588
#> 2       29         T         A         1762       1762            1        588
#> 3      190         C         G         2054       2054            1        685
#> 4      224         A         T         1964       1964            1        655
#> 5       21         G         G         1774       1774            1        592
#> 6       21         G         G         1774       1774            1        592
#>   QUERYID TXID CDSID GENEID   CONSEQUENCE REFCODON VARCODON REFAA VARAA
#> 1     840  116   110   UL54 nonsynonymous      GAC      AAC     D     N
#> 2     840  116   110   UL54 nonsynonymous      GAC      AAC     D     N
#> 3     834  116   110   UL54 nonsynonymous      AAT      AGT     N     S
#> 4     836  116   110   UL54 nonsynonymous      TCA      TTA     S     L
#> 5    1432   56    53   UL97 nonsynonymous      TGC      GGC     C     G
#> 6    1432   56    53   UL97 nonsynonymous      TGC      GGC     C     G
#>   aachange MUTATION_ID VIRUS         GENOTYPE GENE AA_CHANGE Ganciclovir
#> 1    D588N          45  HCMV            AD169 UL54     D588N           2
#> 2    D588N          46  HCMV                  UL54     D588N         3.8
#> 3    N685S         515  HCMV Clinical Isolate UL54     N685S Susceptible
#> 4    S655L         509  HCMV Clinical Isolate UL54     S655L Susceptible
#> 5    C592G         167  HCMV                  UL97     C592G         2.9
#> 6    C592G         430  HCMV                  UL97     C592G           3
#>   Cidofovir Foscarnet Brincidofovir Letermovir Tomeglovir GW275175 Maribavir
#> 1       1.3       2.8            NA                    NA       NA          
#> 2       2.7     3.2-9            NA                    NA       NA          
#> 3                                NA                    NA       NA          
#> 4                                NA                    NA       NA          
#> 5                                NA                    NA       NA          
#> 6                                NA                    NA       NA          
#>   Cyclopropavir                         REF_REVIEW
#> 1            NA                                   
#> 2            NA Rev. Med. Virol. 2016; 26: 161–182
#> 3            NA                                   
#> 4            NA                                   
#> 5            NA Rev. Med. Virol. 2016; 26: 161–182
#> 6           3.3                                   
#>                                                                                               REF_TITLE
#> 1 Phenotypic Diversity of Cytomegalovirus DNA Polymerase Gene Variants Observed after Antiviral Therapy
#> 2                                                                                     [51,67,72,83,122]
#> 3                                                    Antiviral Drug Resistance of Human Cytomegalovirus
#> 4                                                    Antiviral Drug Resistance of Human Cytomegalovirus
#> 5                                                                                                      
#> 6                 Cytomegalovirus UL97 Mutations Affecting Cyclopropavir and Ganciclovir Susceptibility
#>                                                REF_LINK
#> 1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059355/
#> 2                                                      
#> 3 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2952978/
#> 4 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2952978/
#> 5                                                      
#> 6                  https://aac.asm.org/content/55/1/382
#>                                     REF_DOI
#> 1            doi: 10.1016/j.jcv.2011.01.004
#> 2                                          
#> 3                 doi: 10.1128/CMR.00009-10
#> 4                 doi: 10.1128/CMR.00009-10
#> 5                                          
#> 6 https://dx.doi.org/10.1128%2FAAC.01259-10
#>                                                          TEST_METHOD  TM_CLASS
#> 1                                              Recombinant BAC, SEAP  in vitro
#> 2                                    Lab derived mutant, PRA or SEAP  in vitro
#> 3 Observed in baseline sequences or drug-sensetive clinical isolates Anecdotal
#> 4 Observed in baseline sequences or drug-sensetive clinical isolates Anecdotal
#> 5                                    Lab derived mutant, PRA or SEAP  in vitro
#> 6                                              Recombinant BAC, SEAP  in vitro
#>   co.gene co.AA Created_date Created_by note Status
#> 1                 16/03/2020  OJCharles           A
#> 2                 17/11/2019  OJCharles           A
#> 3                 16/03/2020  OJCharles           A
#> 4                 16/03/2020  OJCharles           A
#> 5                 17/11/2019  OJCharles           A
#> 6                 22/11/2019  OJCharles           A
```

#### View all mutations

It may be the case that your sequences contain non synonymous mutations
in resistance genes, that may confer changes in resultant protein.
Although there are a large number of mutations accounted for in this
databse, there is the possiblity to observe something novel. The
call\_resistance function allows for these sorts of questions.

``` r

# returns all variants, with any resistance data, with annotation of protein coded for & if it coded for a change in protein sequence.
mutations_all = call_resistance(infile = my_sample, all_mutations = T)
#> 'select()' returned many:1 mapping between keys and columns


# to view all mutations in resistance genes we can use base R to filter.
mutations_res = mutations_all[mutations_all$GENEID %in% c("UL54", "UL97", "UL27", "UL51", "UL56", "UL89"),]


# are there any non-synonymous (DNA variants that result in a change of amino acid) variants in resistance genes.
mutations_res_nonsyn = mutations_res[mutations_res$CONSEQUENCE == "nonsynonymous",]

# here the top 3 mutations are nonsynonymous, but have no identified resistance effect.
head(mutations_res_nonsyn[,c(1,8,21,32:40)])
#>           change   freq   CONSEQUENCE Ganciclovir Cidofovir Foscarnet
#> 996    UL51_T47M 99.88% nonsynonymous        <NA>      <NA>      <NA>
#> 1002 UL54_A1108T 99.05% nonsynonymous        <NA>      <NA>      <NA>
#> 1004  UL54_A692V   100% nonsynonymous        <NA>      <NA>      <NA>
#> 1007  UL54_D588N   5.2% nonsynonymous           2       1.3       2.8
#> 1008  UL54_D588N   5.2% nonsynonymous         3.8       2.7     3.2-9
#> 1018  UL54_L897S 99.82% nonsynonymous        <NA>      <NA>      <NA>
#>      Brincidofovir Letermovir Tomeglovir GW275175 Maribavir Cyclopropavir
#> 996             NA       <NA>         NA       NA      <NA>            NA
#> 1002            NA       <NA>         NA       NA      <NA>            NA
#> 1004            NA       <NA>         NA       NA      <NA>            NA
#> 1007            NA                    NA       NA                      NA
#> 1008            NA                    NA       NA                      NA
#> 1018            NA       <NA>         NA       NA      <NA>            NA
```

There are a number of variants above in resistance genes, that result in
a change of protein sequence and may deserve attention.

## Run the shiny app

``` r
# runShinyCMV()
```
