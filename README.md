
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cmvdrg

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/ucl-pathgenomics/cmvdrg.svg?branch=master)](https://travis-ci.com/ucl-pathgenomics/cmvdrg)
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
(whole genomes & fragments) which will be mapped to RefSeq NC\_006273.2.
NGS variant data assembled to NC\_006273.2 is accepted in VCF \>= ver4.0
& Varscan2 tab formats.

#### Database

Contains the relationships between:

  - Mutation
  - Drug susceptibility, values are EC50 fold changes to susceptible
    strains. “Susceptible” or “Resistant” are present if data is
    anecdotal.
  - The method for Resistance Phenotyping, including where possible
    strain information used in marker transfer. i.e. AD169, TOWNE.
  - Publication reference, DOI, web address & review paper reference if
    used.
  - Where a co-mutation susceptibility profile is available for a known
    drug, this is included.
  - Administrative information, last update to row & whether it is used
    in calculations or not. i.e. A “Active”" is included, where there is
    ambiguity the status field is “U” Unused, or “R” to Review

#### Web service

A user-friendly Shiny Applications has been bundled with this package.
The same application is available over the internet here
<http://cmv-resistance.ucl.ac.uk/cmvdrg/> where the terms of use are
contained.

## Installation

You can install the current version from
[GitHub](https://github.com/ucl-pathgenomics/cmvdrg) with:

``` r
# install.packages("devtools")
devtools::install_github("ucl-pathgenomics/cmvdrg")
```

Dependencies for FASTA file handling are MAFFT and SNP-Sites available
preferably via conda. snp-sites \>= 2.3 has been tested.

``` bash
conda config --add channels bioconda
conda install snp-sites
conda install mafft
```

## Usage

``` r
library("cmvdrg")

## call resistant variants


my_sample = system.file("testdata", "A10.vcf", package = "cmvdrg")


data = call_resistance(infile = my_sample, all_mutations = F)


print(data[,c("change", "freq", "Ganciclovir")])
#> [1] change      freq        Ganciclovir
#> <0 rows> (or 0-length row.names)




## call all variants

mutations_all = call_resistance(infile = my_sample, all_mutations = T)


#to view all mutations in resistance genes we can filter
mutations_res = mutations_all[mutations_all$GENEID %in% c("UL54", "UL97", "UL27", "UL51", "UL56", "UL89"),]


# are there any non-synonymous (DNA variants that result in a change of amino acid) variants in resistance genes
mutations_res_nonsyn = mutations_res[mutations_res$CONSEQUENCE == "nonsynonymous",]


# here the top 3 mutations are nonsynonymous, with no identified resistance effect.
head(mutations_res_nonsyn[,c(1,8,21,32:40)])
#>          change   freq   CONSEQUENCE Ganciclovir Cidofovir Foscarnet
#> 996   UL51_A47V 99.88% nonsynonymous        <NA>      <NA>      <NA>
#> 997  UL51_Q112H   100% nonsynonymous        <NA>      <NA>      <NA>
#> 1004 UL54_A588T   5.2% nonsynonymous        <NA>      <NA>      <NA>
#> 1005 UL54_A692V   100% nonsynonymous        <NA>      <NA>      <NA>
#> 1006 UL54_A695G  1.05% nonsynonymous        <NA>      <NA>      <NA>
#> 1009 UL54_E939D 99.78% nonsynonymous        <NA>      <NA>      <NA>
#>      Brincidofovir Letermovir Tomeglovir GW275175 Maribavir Cyclopropavir
#> 996             NA       <NA>         NA       NA      <NA>            NA
#> 997             NA       <NA>         NA       NA      <NA>            NA
#> 1004            NA       <NA>         NA       NA      <NA>            NA
#> 1005            NA       <NA>         NA       NA      <NA>            NA
#> 1006            NA       <NA>         NA       NA      <NA>            NA
#> 1009            NA       <NA>         NA       NA      <NA>            NA


## run the shiny application
# runShinyCMV()
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on the GitHub Issues page. For questions and other
discussions feel free to contact. [Oscar Charles -
maintainer](mailto:oscar.charles.18@ucl.ac.uk)
