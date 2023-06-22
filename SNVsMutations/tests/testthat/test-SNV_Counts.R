### Load the necessary packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(VariantAnnotation)

### Create the necessary variables to test the functions
fl <- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")
vcf <- readVcf(fl, "hg38")[1:10]
ref_genome <- Hsapiens
mut_type <- SNV_Types(vcf, context_length = 5, ref_genome = ref_genome)

### Input check
test_that("SNV_Counts throws an error if mutationstype_vector is not the correct format", {
  a <- c("[A>T]", "C>G", "G>T")
  expect_error(SNV_Counts(mutationstype_vector = a))
})

test_that("SNV_Counts throws an error if both mutationstype_vector and vcf files are provided", {
  expect_error(SNV_Counts(mutationstype_vector = mut_type, vcf_file=vcf))
})


test_that("SNV_Counts throws an error if required arguments are missing when mutationstype_vector is not provided", {
  expect_error(SNV_Counts(context_length=7))
})


### Output check
test_that("SNV_Counts computes mutation types correctly from mutationstype_vector", {
  count_table <- SNV_Counts(mutationstype_vector = mut_type, onlySNV=TRUE)
  expect_equal(nrow(count_table), 2)
  expect_equal(as.vector(count_table$Mutation), c("C>T", "T>C"))
})

test_that("SNV_Counts computes mutation types correctly from VCF file, reference genome, and context length", {
  count_table <- SNV_Counts(vcf_file = vcf, ref_genome = ref_genome, context_length = 7, onlySNV=TRUE)
  expect_true(all(grepl("[.*>.*]", count_table$Mutation)))
})
