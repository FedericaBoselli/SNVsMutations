# SNVsMutations
Bioconductor compliant R package that offers functions for identifying the category of single nucleotide variants (SNVs) contained within VCF files.


The package analyzes a set of single nucleotide variants (SNVs) in a VCF file, along with the associated reference genome, and utilizes a parameter called "context_length" to identify and classify the mutation type for each individual mutation. The package provides flexibility by allowing users to either provide a mutation file or calculate mutation types directly from the provided data. It is also able to provide a graphical representation of the mutation counts for each specific mutation type. 
