#' Get fasta sequences for a gene list.
#'
#' @description Get fasta sequences of given input list from reference
#'   multi-fasta file.
#' @param gene_list A charcter vector of gene names or gene ids;
#' the gene names or gene ids should exactly match with sequence
#' headers present in multi-fasta file
#' @param fasta_file Either a path or a connection to reference multi-fasta
#'   file, from which subset of sequences for given input list is to be
#'   retrieved.
#' @param outfile A character vector containing the path to the file to write
#'   output
#' @importFrom Biostrings readBStringSet writeXStringSet
#' @importFrom tidyr tibble
#' @importFrom dplyr mutate
#' @import magrittr
#' @return A multi-fasta file, containing sequences of given input list id's.
#' @export
#'
#' @examples
#' \dontrun{
#'
#'  myGenelist <- system.file("exdata", "Sc_myGenelist.txt", package = "fastaR")
#'  myGenelist <- scan(myGenelist,  what="character", sep=NULL)
#'
#'  ref_fasta <- system.file("exdata", "Sc_nucl_R64-2-1.fasta", package = "fastaR")
#'  fastaR::fa_some_records(gene_list=myGenelist, fasta_file=ref_fasta, outfile="sc_myGenelist.fa")
#'
#' }
fa_some_records <- function(gene_list, fasta_file, outfile="stdout.fa"){

          ##Read fasta file
          input_sequence = Biostrings::readBStringSet(filepath=fasta_file,use.names = TRUE)
         # names(input_sequence) <- base::gsub('[[:space:]| :].*', '', names(input_sequence))

          renamed_seq <- names(input_sequence) %>%
                              tidyr::tibble(rownames = .) %>%
                              dplyr::mutate(seq_new=base::gsub('[[:space:]| :].*', '', rownames))

          ID <- as.matrix(gene_list)
          overlapped_id <- subset(renamed_seq, renamed_seq$seq_new %in% ID)

          #message("no. of ID's: ", nrow(ID))

          ## Find overlap between ID's and fasta file

          overlap <-  input_sequence[which(names(input_sequence) %in% overlapped_id$rownames),]

          ## write fasta file

          Biostrings::writeXStringSet(x =overlap, filepath = outfile, format = "fasta")

          return(overlap)

}


#' Get size of fasta sequence
#' @description Get size (in bp) for each sequence in multi-fasta file.
#'
#' @param fasta_file Either a path or a connection to multi-fasta file.
#'   The input sequence file should have extention .fa or .fasta
#'   In the sequence header: only string before first space and/or first colon (:) will be considered for futher processes.
#'   **Important consideration when header have big names.
#'
#' @return A tibble of gene_id and gene length
#' @export
#' @importFrom Biostrings readBStringSet width
#' @importFrom dplyr bind_cols
#' @examples
#' \dontrun{
#'
#'  ref_fasta <- system.file("exdata", "Sc_nucl_R64-2-1.fasta", package = "fastaR")
#'
#'  fastaR::fa_size(fasta_file=ref_fasta)
#'
#' }
fa_size <- function(fasta_file){

          #---- read input genome fasta sequence
          input_sequence <- Biostrings::readBStringSet(filepath=fasta_file,use.names = TRUE)

          #---- extract names and chromosome length
          names(input_sequence) <- base::gsub('[[:space:]| :].*', '', names(input_sequence))

          genome_size_mat <- dplyr::bind_cols(Seq_id=base::names(input_sequence),Length= Biostrings::width(input_sequence))

          #colnames(genome_size_mat) <- c("Seq_id", "Length")

          #message( paste(output_name,".size", sep=""), " file is saved in working directory!")

          return(genome_size_mat)
}

#' Fasta file summary
#' @description Get summary of input multi-fasta file like mean length, median
#'   length, GC content etc
#' @param fasta_file Either a path or a connection to multi-fasta file.
#'   The input sequence file should have extention .fa or .fasta
#'   In the sequence header: only string before first space and/or first colon (:) will be considered for futher processes.
#'   **Important consideration when header have big names.
#' @return A summary table
#' @export
#' @importFrom Biostrings readBStringSet letterFrequency width
#' @importFrom dplyr bind_cols rename summarise n mutate
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling row_spec column_spec
#'
#' @examples
#' \dontrun{
#'
#'  ref_fasta <- system.file("exdata", "Sc_nucl_R64-2-1.fasta", package = "fastaR")
#'  fastaR::fa_summary(fasta_file=ref_fasta)
#' }
fa_summary <- function(fasta_file){

                    # read input sequences as biostring object
                    input_sequence <- Biostrings::readBStringSet(filepath=fasta_file,use.names = TRUE)
                    base::names(input_sequence) <- base::gsub('[[:space:]| :].*', '', base::names(input_sequence))

                    # compute input sequence lengths
                    total_length <-  Biostrings::letterFrequency(input_sequence,letters = "ATGC", collapse=TRUE)
                    gc_content <-  Biostrings::letterFrequency(input_sequence, "GC", collapse=TRUE)


                    summary_length <- dplyr::bind_cols(Seq_id=base::names(input_sequence), Length=Biostrings::width(input_sequence)) %>%
                                        dplyr::summarise( num_of_seq=dplyr::n(),min=min(Length), max=max(Length),mean=round(mean(Length),2),
                                                   median= round(median(Length),2)) %>%
                                        dplyr::mutate(percent_gc = round(100*(gc_content/total_length),2))


                    knitr::kable(t(summary_length), col.names = c("Summary"),"html", align = "l") %>%
                              kableExtra::kable_styling(bootstrap_options = c("striped", "condensed", "responsive"), full_width = F,font_size =14,  stripe_color = "aquamarine3") %>%
                              kableExtra::row_spec(0,bold = TRUE, italic = FALSE, color = "black") %>%
                              kableExtra::column_spec(1:2, bold=FALSE, color="blue")

                    return(summary_length)


}


#' Get GC percent for multi-fasta file
#' @description for given multi-fasta file get GC percent for each sequence
#' @param fasta_file Either a path or a connection to multi-fasta file. The
#'   input sequence file should have extention .fa or .fasta
#'   In the sequence header: only string before first space and/or first colon (:) will be considered for futher processes.
#'   **Important consideration when header have big names.
#' @return A tibble of GC percent and a barplot of GC percent distribution
#'  \code{if the fasta file contains sequences less than 20}
#' @export
#' @importFrom Biostrings readBStringSet letterFrequency
#' @importFrom tidyr as_tibble
#' @importFrom dplyr mutate select
#' @importFrom ggplot2 ggplot aes geom_col coord_flip theme_classic geom_text theme element_blank element_text
#' @examples
#' \dontrun{
#'  ref_fasta <- system.file("exdata", "Sc_nucl_R64-2-1.fasta", package = "fastaR")
#'  fastaR::fa_percent_GC(fasta_file=ref_fasta)
#' }
fa_percent_GC <- function(fasta_file){

          # read input sequences as biostring object
          input_sequence <- Biostrings::readBStringSet(filepath=fasta_file,use.names = TRUE)
          base::names(input_sequence) <- base::gsub('[[:space:]| :].*', '', base::names(input_sequence))

          output_name <- base::gsub(".fa.*","", basename(fasta_file))

          # compute input sequence lengths
          total_length <-  Biostrings::letterFrequency(input_sequence,letters = "ATGC", collapse=FALSE)
          gc_content <-  Biostrings::letterFrequency(input_sequence, "GC", collapse=FALSE)

          nucl_table <- cbind(names=base::names(input_sequence), total_length, gc_content) %>% tidyr::as_tibble()
          base::colnames(nucl_table) <- c("names","sequence_len", "gc")

          gc_each_seq <- nucl_table %>%
                              dplyr::mutate(percent_gc = round(100*(as.numeric(gc)/as.numeric(sequence_len)),2)) %>%
                              dplyr::select(names, percent_gc)

          if(nrow(gc_each_seq) <= 20){

                    message("plotting GC percent for each sequence")

                    gg <- gc_each_seq %>%
                              ggplot2::ggplot(ggplot2::aes(names, percent_gc, fill=names, label=percent_gc))+
                              ggplot2::geom_col(alpha=0.8)+
                              ggplot2::coord_flip()+
                              ggplot2::theme_classic()+
                              ggplot2::geom_text(hjust=1, size=5)+
                              ggplot2::theme(axis.title.y=ggplot2::element_blank(),
                                             axis.text=ggplot2::element_text(color = "black", size = 12),
                                             legend.position = "none")
                    print(gg)
                    return(gc_each_seq)
          }

          else{
                    return(gc_each_seq)
          }

}
