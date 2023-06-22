#' SNV_Types Function
#'
#' This function extracts SNV (Single Nucleotide Variant) mutations from a VCF 
#' file and generates mutation strings with reference and alternate alleles. 
#' It also considers a context length to extract the upstream and downstream bases flanking the mutation.
#'
#' @param vcf The loaded VCF (Variant Call Format) file containing the mutations.
#' @param context_length The length of the context (number of bases) to consider 
#' around (upstream and downstream) each mutation. It must be an odd positive number.
#' @param ref_genome The corresponding reference genome (e.g., human genome hg38).
#'
#' @return A vector containing each SNV mutation in the VCF file.
#' @import S4Vectors
#' @import IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges seqnames
#' @importFrom MatrixGenerics rowRanges
#' @importFrom Biostrings reverseComplement
#' @import VariantAnnotation
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @examples 
#' f <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' vcf <- readVcf(f, "hg38")[1:100]
#' ref_genome <- Hsapiens
#' SNV_Types(vcf, context_length = 3, ref_genome)   
#' 
#' @export

SNV_Types <- function(vcf, context_length, ref_genome) {
  
  ## context length must always be odd!
  if (context_length %% 2 == 0) {
    stop("context_length must be odd")
  }
  
  # filter VCF only to SNV mutations
  vcf_snv <- subset(vcf, isSNV(vcf))
  
  # retrieve the genomic positions of single-nucleotide variants from the vcf_snv file
  snv <- rowRanges(vcf_snv)
  
  # calculate the number of single-nucleotide variants stored in snv
  num_snv <- length(snv)
  
  # setting the up_down variable for the next operations
  up_down <- (context_length - 1) / 2 
  
  # create empty vector to store the mutations
  mutations <- vapply(seq_len(num_snv), function(i) {
    chr <- as.character(seqnames(snv[i]))
    snv_start <- start(snv[i])
    
    # set start and end of the sequence according to context_length parameter
    start <- (snv_start - up_down)
    end <- (snv_start + up_down)
    
    # create a genomic range object with chromosome information and start-end positions
    gr <- GRanges(seqnames = paste0("chr", chr), ranges = IRanges(start = start, end = end))
    
    # retrieve the corresponding reference sequence from the provided reference genome
    reference_seq <- getSeq(ref_genome, gr)
    
    # extract the reference allele and alternate alleles from the snv object at index i
    ref <- snv[i]$REF
    alt <- unlist(snv[i]$ALT)
    
    # Check if the reference allele (ref) is either "C" or "T". 
    # If it is, the nucleotides, ref, and alt variables are assigned the 
    # corresponding values from reference_seq, ref, and alt.
    if (as.character(ref) %in% c("C", "T")) {
      nucleotides <- as.character(reference_seq)
      ref <- as.character(ref)
      alt <- as.character(alt)
    } else {
      # Otherwise, the code assigns the reverse complement values of 
      # reference_seq, ref, and alt to nucleotides, ref, and alt, respectively.
      nucleotides <- as.character(reverseComplement(reference_seq))
      ref <- as.character(reverseComplement(ref))
      alt <- as.character(reverseComplement(alt))
    }
    
    # split the reference sequence in UPSTREAM and DOWNSTREAM bases
    UP_bases <- substr(reference_seq, 1, up_down)
    DOWN_bases <- substr(reference_seq, (context_length - up_down) +1, context_length)
    
    # Construct mutation string
    mutation <- paste0(UP_bases, "[", ref, ">", alt, "]", DOWN_bases) 
    
    # Return mutation string
    mutation
  }, character(1))
  
  # Return vector of mutations
  return(mutations)
}










