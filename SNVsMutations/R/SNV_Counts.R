#' SNV_Counts Function
#'
#' This function counts the occurrences of mutation types and creates a 
#' visualization plot based on the provided mutation file or by computing the 
#' mutation types from a VCF file, reference genome, and context length.
#'
#' @param mutationstype_vector A vector of mutation types obtained from SNV_Types function
#'        or having the same structure.
#' @param vcf_file the VCF file.
#' @param ref_genome Reference genome object.
#' @param context_length Length of the context (number of bases) to consider 
#' around (upstream and downstream) each mutation. It must be an odd positive number. 
#' @param plot Logical value indicating whether to create a visualization plot.
#' @param onlySNV Logical value indicating whether to count only SNV mutations.
#' @return A data frame containing the counts of each SNV mutation type.
#' @details If a mutationstype_vector is provided, it is used directly to count the mutation types.
#' Otherwise, the mutation types are computed using the SNV_Counts function 
#' based on the vcf_file, ref_genome, and context_length parameters. 
#' The function also creates a bar plot showing the mutation counts with 
#' different colors for each mutation type.
#' The plot is displayed if plot is set to TRUE. In particular:
#' If onlySNV is TRUE, the plot shows the mutation types by considering only 
#' reference and alternative allele (SNV only).
#' If onlySNV is FALSE, the plot shows the mutation types by also considering 
#' the bases upstream and downstream (taking into account the 'context_length' variable) each SNV. 
#' 
#' @import ggplot2
#' @importFrom grDevices hcl
#' @importFrom stats setNames
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 scale_fill_identity
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 ggtitle
#' 
#' @examples
#' f <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' vcf <- readVcf(f, "hg38")[1:100]
#' ref_genome <- Hsapiens
#' mut_types <- SNV_Types(vcf, context_length = 3, ref_genome)
#' SNV_Counts(mut_types, plot=TRUE, onlySNV=TRUE)   
#' 
#' @export


SNV_Counts <- function(mutationstype_vector = NULL, vcf_file = NULL, 
                       ref_genome = NULL, context_length = NULL,
                       plot = TRUE, onlySNV = FALSE){
  
  # If the file of mutation types is provided:
  if(!is.null(mutationstype_vector)) {
    if(!(length(mutationstype_vector) == sum(grepl('\\[[CT]>[AGCT]\\]', 
                                                   mutationstype_vector)))) {
      stop('Error: the mutation file must be obtained via SNV_Types function 
           (or must have the same structure)!')
    }
    if(any(!is.null(vcf_file), !is.null(ref_genome), !is.null(context_length))) {
      stop('If you provide the mutation file, you can disregard the VCF file, 
           reference genome, and context length inputs.')
    }
    if(!onlySNV){
      Mutation <- mutationstype_vector
    } else {
      # The code uses the gsub function to extract the mutation types from 
      # mutationstype_vector by removing all text except the content within 
      # square brackets and assigning the results to SNVs_types. The gsub function 
      # performs a global search and replaces operation based on a specified pattern.
      Mutation <- gsub(".*\\[(.*?)\\].*", "\\1", mutationstype_vector)
    }
  }
  
  # If there is no mutations type file, the types are computed with the function
  # SNV_Types from a vcf file, the reference genome and the context length parameter:
  else {
    if(any(is.null(vcf_file), is.null(ref_genome), is.null(context_length))){
      stop('You must provide either a mutation file or include the reference genome, 
           VCF file, and context length together.')}
    mutationstype_vector <- SNV_Types(vcf_file, context_length, ref_genome)
    if(!onlySNV){
      Mutation <- mutationstype_vector
    } else {
      Mutation <- gsub(".*\\[(.*?)\\].*", "\\1", mutationstype_vector)
    }
  }
  
  # CREATE BAR PLOT FOR VISUALIZATION
  count_table <- as.data.frame(table(Mutation))
  
  if (plot) {
    # Extract values between squared brackets and generate pastel colors 
    # based on unique values inside the brackets
    mutation_values <- gsub(".*\\[(.*?)\\].*", "\\1", count_table$Mutation)
    unique_values <- unique(mutation_values)
    num_colors <- length(unique_values)
    colors <- hcl(h = seq(15, 375, length = num_colors + 1), c = 30, l = 70)
    color_map <- setNames(colors, unique_values)
    Color <- color_map[mutation_values]
    count_table <- cbind(count_table, Color)
    
    # Create the plot and assign a unique color to each bar
    g <- ggplot(count_table, aes(x = Mutation, y = count_table$Freq, 
                                 fill = Color)) + 
      geom_bar(stat = 'identity') + 
      coord_flip() +
      labs(y="Frequency", title ="Mutation Counts") +
      scale_fill_identity() +
      theme_minimal()
    print(g)
  }
  return(count_table)
}

