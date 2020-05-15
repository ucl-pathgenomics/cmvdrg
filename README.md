# cmvdrg - An R package for Human Cytomegalovirus antiviral Drug Resistance Genotyping


## R package

 1. Install devtools if needed with: install.packages(“devtools”)
 2. Install cmvdrg by running: devtools::install_github("ucl-pathgenomics/cmvdrg”). This command will take care of installing other dependencies. You may need to upgrade your version of R.
 3. A vignette has been provided which acts as a tutorial. Test sequences have been provided.

## Database
A manually curated csv file, containing the relationships between 

 - Gene
 - AA mutation
 - Drug susceptibility, values are EC50 fold changes to susceptible strains. "Susceptible" or "Resistant" are present if data is anecdotal.
 - The method for Resistance Phenotyping, including where possible strian information used in marker transfer. i.e. AD169, TOWNE.
 - Publication reference, DOI, web address & review paper reference if used.
 - Where a multi-gene, multi AA mutation susceptibility profile is available for a known drug, this is included too.
 - Administrative information, last update to reference, updater, whether the row is used in calculations or not. i.e. if a mutations is present in a review paper, but is not traceable to given references this column is set to "U" Unused, or "R" to Review.

## Web service
A user friendly Shiny Applications has been bundled with this package. The same application is available over the internet here [http://51.11.13.133:3838/cmvdrg/](http://51.11.13.133:3838/cmvdrg/) where the terms of use are also contained. 
