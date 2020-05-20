context("calling resistance")

library(cmvdrg)


test_that("variant files return resistance", {
  expect_equal(nrow(call_resistance(infile = system.file("testdata",  "A10.vcf", package = "cmvdrg"))), 11)
  expect_equal(nrow(call_resistance(infile = system.file("testdata",  "A10.vcf", package = "cmvdrg"),all_mutations = T)), 1980)
  expect_equal(nrow(call_resistance(infile = system.file("testdata",  "A10.tab", package = "cmvdrg"))), 4)
  expect_equal(nrow(call_resistance(infile = system.file("testdata",  "snpsites.vcf", package = "cmvdrg"))), 2)
})

# travis seems a reall faff to get working for snp-sites 2.3 as it's r builds are xenial not bionic ubuntu. so ignoring fasta tests for now.
# tested regularly locally. mafft and snp-sites don't need to be tested alone.



