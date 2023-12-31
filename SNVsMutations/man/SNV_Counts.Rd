% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SNV_Counts.R
\name{SNV_Counts}
\alias{SNV_Counts}
\title{SNV_Counts Function}
\usage{
SNV_Counts(
  mutationstype_vector = NULL,
  vcf_file = NULL,
  ref_genome = NULL,
  context_length = NULL,
  plot = TRUE,
  onlySNV = FALSE
)
}
\arguments{
\item{mutationstype_vector}{A vector of mutation types obtained from SNV_Types function
or having the same structure.}

\item{vcf_file}{the VCF file.}

\item{ref_genome}{Reference genome object.}

\item{context_length}{Length of the context (number of bases) to consider
around (upstream and downstream) each mutation. It must be an odd positive number.}

\item{plot}{Logical value indicating whether to create a visualization plot.}

\item{onlySNV}{Logical value indicating whether to count only SNV mutations.}
}
\value{
A data frame containing the counts of each SNV mutation type.
}
\description{
This function counts the occurrences of mutation types and creates a
visualization plot based on the provided mutation file or by computing the
mutation types from a VCF file, reference genome, and context length.
}
\details{
If a mutationstype_vector is provided, it is used directly to count the mutation types.
Otherwise, the mutation types are computed using the SNV_Counts function
based on the vcf_file, ref_genome, and context_length parameters.
The function also creates a bar plot showing the mutation counts with
different colors for each mutation type.
The plot is displayed if plot is set to TRUE. In particular:
If onlySNV is TRUE, the plot shows the mutation types by considering only
reference and alternative allele (SNV only).
If onlySNV is FALSE, the plot shows the mutation types by also considering
the bases upstream and downstream (taking into account the 'context_length' variable) each SNV.
}
\examples{
f <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(f, "hg38")[1:100]
ref_genome <- Hsapiens
mut_types <- SNV_Types(vcf, context_length = 3, ref_genome)
SNV_Counts(mut_types, plot=TRUE, onlySNV=TRUE)   

}
