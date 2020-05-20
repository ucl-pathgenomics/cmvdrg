context("calling resistance")

library(cmvdrg)


test_that("variant files return resistance", {
  expect_equal(nrow(call_resistance(infile = system.file("testdata",  "A10.vcf", package = "cmvdrg"))), 11)
  expect_equal(nrow(call_resistance(infile = system.file("testdata",  "A10.vcf", package = "cmvdrg"),all_mutations = T)), 1980)
  expect_equal(nrow(call_resistance(infile = system.file("testdata",  "A10.tab", package = "cmvdrg"))), 4)
  expect_equal(nrow(call_resistance(infile = system.file("testdata",  "snpsites.vcf", package = "cmvdrg"))), 2)
})



test_that("fasta files return resistance", {
  expect_equal(nrow(call_resistance(infile = system.file("testdata",  "A10.fasta", package = "cmvdrg"))), 2)
  expect_equal(nrow(call_resistance(infile = system.file("testdata",  "A10_fragment.fasta", package = "cmvdrg"))), 2)
})
