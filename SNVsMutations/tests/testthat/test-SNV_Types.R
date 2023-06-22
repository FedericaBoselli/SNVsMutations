### Load the necessary packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(VariantAnnotation)

### Input check 

test_that("SNV_Types throws an error for even context_length", {
  fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  vcf <- readVcf(fl, "hg38")[1:20]
  context_length <- 4  # even context_length
  ref_genome <- Hsapiens
    
    expect_error(SNV_Types(vcf, context_length, ref_genome))
})


### Output check

test_that("SNV_Types filters only SNV mutations from the VCF file", {
  fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  vcf <- readVcf(fl, "hg38")[1:20]
  context_length <- 5  # even context_length
  ref_genome <- Hsapiens
    
    mutations <- SNV_Types(vcf, context_length, ref_genome)
  
  expect_true(all(grepl("\\[[ACTG]>[ACTG]\\]", mutations)))
})


test_that("SNV_Types returns the correct format for the mutations with a given context length", {
  fl <- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")
  vcf <- readVcf(fl, "hg38")[1:20]
  context_length <- 5
  ref_genome <- Hsapiens
  
  mutations <- SNV_Types(vcf, context_length, ref_genome)
  
  # Check if all mutations have the correct format
  for (mutation in mutations) {
    expect_true(grepl("^..\\[[ACGT]>[ACGT]\\]..$", mutation))
  }
})
