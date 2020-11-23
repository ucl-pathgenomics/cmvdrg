#' Returns a dataframe of the cmvdrg database
#' @return dataframe of the cmvdrg database
#' @export

cmvdrg_data <- function(){
  
  file = system.file("db", "cmvdrg-db1.csv", package = "cmvdrg")
  dat = utils::read.csv(file, stringsAsFactors = F)
  colnames(dat)[[1]] = "mutation_id"
  return(dat)
}
