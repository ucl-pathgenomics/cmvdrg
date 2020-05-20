
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
<!-- badges: end -->

## overview

cmvdrg is a R package to enable antiviral drug resistance genotyping,
with Human Cytomegalovirus sequencing data. Accepted inputs are FASTA
whole genomes, FASTA fragments which will be mapped to Refseq
NC\_006273.2. NGS variant data assembled to NC\_006273.2 is accepted in
VCF \>= ver4.0 & Varscan2 tab formats.

#### Database

Contains the relationships between:

  - Mutation
  - Drug susceptibility, values are EC50 fold changes to susceptible
    strains. “Susceptible” or “Resistant” are present if data is
    anecdotal.
  - The method for Resistance Phenotyping, including where possible
    strian information used in marker transfer. i.e. AD169, TOWNE.
  - Publication reference, DOI, web address & review paper reference if
    used.
  - Where a co-mutation susceptibility profile is available for a known
    drug, this is included.
  - Administrative information, last update to reference, updater,
    whether the row is used in calculations or not. i.e. if a mutations
    is present in a review paper, but is not traceable to given
    references this column is set to “U” Unused, or “R” to Review.

#### Web service

A user friendly Shiny Applications has been bundled with this package.
The same application is available over the internet here
<http://51.11.13.133:3838/cmvdrg/> where the terms of use are contained.

## Installation

You can install the current version from
[GitHub](https://github.com/ucl-pathgenomics/cmvdrg) with:

``` r
# install.packages("devtools")
devtools::install_github("ucl-pathgenomics/cmvdrg")
```

Dependendencies for FASTA file handling are MAFFT and SNP-Sites
available preffereably via conda. snp-sites \>= 2.3 has been tested.

``` bash
conda config --add channels bioconda
conda install snp-sites
conda install mafft
```

## Usage

``` r
library("cmvdrg")

## call resistance mutations only
my_sample = system.file("testdata", "A10.vcf", package = "cmvdrg")


# this calls any variants in the sample file, that are also present in the cmvdrg database.
data = call_resistance(infile = my_sample, all_mutations = F)


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



## call all mutations - i.e. there may be interesting non-synonymous mutants predsent in resistance genes.
# returns all variants, with any resistance data, with annotation of protein coded for & if it coded for a change in protein sequence.
mutations_all = call_resistance(infile = my_sample, all_mutations = T)


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


## run the shiny application
# runShinyCMV()
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on the GitHub Issues page. For questions and other
discussions feel free to contact. [Oscar Charles -
maintainer](mailto:oscar.charles.18@ucl.ac.uk)
