#' Clinical information table
#'
#' Produces a data.table with a clinical overview for resistance phenotype against all drugs in that sample
#'
#' @param f.dat the complete resistance data.frame
#' @return a data.table
#' @export

make_clin_table = function(f.dat){

  # this function processes the aligned genotypic & resistance data and creates outputs.
  # as a result there are a few decisions to be made here, where the code acts as an arbitrator (not good)
  # this is overcome, by also providing all the raw data in another output.
  
  # data wrangling
  # manually set
  dat_drug = f.dat[,c(32:40)]
  
  
  # set up the empty dataframe 
  #dat = data.frame(CDV = c(4,"strong"), FOS  = c(1,"strong"),  GCV = c(1,"anecdotal"), Maribavir =c(0,"none"))
  dat = data.frame(matrix(nrow = 2, ncol = ncol(dat_drug)))
  colnames(dat) = colnames(dat_drug)
  rownames(dat) = c("Resistance Phenotype", "Evidence Strength")
  dat_drug = f.dat[,c(32:40, which(names(f.dat)=="TM_CLASS"))]
  
  # write the table
  # for each drug, see if there are numbers (fold changes are best quality, then in vitro res/sus then anything else)
  for(col in 1:ncol(dat)){
  #for(col in 1:1){ 
   res.pheno = "NA"
    res.ev = "No Evidence"
    col.name = colnames(dat[col])
    t.dat = dat_drug[,c(col, ncol(dat_drug))]
    
    # fix any reference data points where there is a numeric range of fold change ratio values. take lowest value - again arbitration
    if(length(t.dat[base::grepl(pattern = "-",x = t.dat[,1]),1]) > 0){
      t.dat[base::grepl(pattern = "-",x = t.dat[,1]),1] = stringr::str_split(t.dat[base::grepl(pattern = "-",x = t.dat[,1]),1], "-", simplify = T)[,1]
    }
    # are there numbers
    t.grep = base::grepl(pattern = "[0-9]", t.dat[,1])
    if(length(t.grep[t.grep==TRUE]) > 0){
      #see numbers
      #categorised as of ---  The Third International Consensus Guidelines on the Management of Cytomegalovirus in Solid-organ Transplantation
      res.pheno = as.numeric(t.dat[base::grepl(pattern = "[0-9]", t.dat[,1]),1])
      res.pheno = max(res.pheno)
      if(res.pheno >= 15){
        res.pheno = "High level"
      }else if(res.pheno <15 && res.pheno >= 5){
        res.pheno = "Moderate level"
      }else if(res.pheno <5 && res.pheno >= 2){
        res.pheno = "Low level"
      }else if(res.pheno < 2){
        res.pheno = "Susceptible"
      }
    res.ev = "good, in vitro"
    #else if there are only 
    }else if(length(base::grepl(pattern = "[a-z]", t.dat[,1])[base::grepl(pattern = "[a-z]", t.dat[,1])==TRUE]) > 0){
      res.pheno = t.dat[base::grepl(pattern = "[a-z]", t.dat[,1]),]
      count_sus = base::grepl(pattern = "usc", res.pheno[,1])
      count_sus = length(count_sus[count_sus==TRUE])
      count_res = base::grepl(pattern = "esis", res.pheno[,1])
      count_res = length(count_res[count_res==TRUE])
      
      if(count_sus > 0 && count_res == 0){
        res.pheno = "Susceptible"
      }else if(count_sus == 0 && count_res > 0){
        res.pheno = "Resistant, magnitude unknown "
      }else{
        res.pheno = "No Consensus"
      }
      
      res.ev = "weaker, anecdotal"
      
    }else{
      #default is NA so do nothing
    }
    
    
    #write Resistance Phenotype
    dat[1,col] = res.pheno
    # Write Evidence Strength
    dat[2,col] = res.ev
    
  }
  
  # now make output
  
  #js <- "(/High/).test(value) ? 'red' : (/Low/).test(value) ? 'yellow' : (/Susc/).test(value) ? 'green' : ''"
  #js <- "(/High/).test(value) ? '#ff6f69' : (/Moderate/).test(value) ? '#ff6f69' : (/Low/).test(value) ? '#ffcc5c' : (/Susc/).test(value) ? '#96ceb4' : (/vitro/).test(value) ? '#96ceb4' : (/anecdotal/).test(value) ? '#ffcc5c' :''"
  js <- "(/High/).test(value) ? '#C34318' : (/Moderate/).test(value) ? '#F68C1B' : (/Low/).test(value) ? '#FFC605' : (/Susc/).test(value) ? '#759F2F' : (/vitro/).test(value) ? '#759F2F' : (/anecdotal/).test(value) ? '#F68C1B' :''"
  
  
  out = DT::datatable(dat, options = list(dom = 't')) %>% 
    DT::formatStyle(names(dat),
                1:ncol(dat), backgroundColor = DT::JS(js))
  
  

  
  
  # out <- datatable(dat, options = list(dom = 't')) %>%
  #   formatStyle(
  #     names(dat),
  #     backgroundColor = styleInterval(c(1, 2), c('white', 'blue', 'red')),
  #     fontWeight = 'bold')

  
  return(out)
}




