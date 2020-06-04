context("displaying visuals")

library(cmvdrg)
library(stringr)

cls <- c("datatables", "htmlwidget")

test_that("Data table is functioning", {
  dat = call_resistance(infile = system.file("testdata",  "A10.vcf", package = "cmvdrg") ,all_mutations = FALSE, inc_anecdotal = T)
  clin_table = make_clin_table(dat)  
  expect_equal(sum(str_count(as.character(clin_table$x$data[1,]),pattern = "Low level")),3)
  expect_equal(sum(str_count(as.character(clin_table$x$data[2,]),pattern = "good, in vitro")),5)
  expect_equal(length(attr(clin_table$x,"colnames")),10)
})


test_that("lollipops are correct", {
  #package variables
  global = list()
  global$res_table = system.file("db", "cmvdrg-db1.csv", package = "cmvdrg")
  #create unique session folder
  global$date <- format(Sys.time(), "%Y-%m-%d")
  global$dir = ""
  
  dat = call_resistance(infile = system.file("testdata",  "A10.vcf", package = "cmvdrg"),all_mutations = FALSE, inc_anecdotal = T)
  plot = plot_lollipop(dat, f.gene = "UL54", global = global)
  
  expect_equal(length(plot$layers),5)
  expect_equal(class(plot$layers[[1]]$geom)[1], "GeomSegment")
  expect_equal(plot$labels$y, "Mutation Frequency")
  
})
