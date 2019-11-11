#' Get fasta sequences for a gene list.
#'
#' @description Get fasta sequences of given input list from reference
#'   multi-fasta file.
#' @param gene_list A charcater vector of gene names; present in multi-fasta
#'   file
#' @param fasta_file Either a path or a connection to reference multi-fasta
#'   file, from which subset of sequences for given input list is to be
#'   retrieved.
#' @param outfile A character vector containing the path to the file to write
#'   output
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings writeXStringSet
#' @return A multi-fasta file, containing sequences of given input list
#' @export
#'
#' @examples
#' \dontrun{
#'
#'  myGenelist <- scan("Sc_myGenelist.txt",  what="character", sep=NULL)
#'
#'  fa_some_records(gene_list=myGenelist, fasta_file="Sc_orf_trans_all_R64-2-1.fasta", outfile="sc_myGenelist.fa")
#'
#' }
fa_some_records <- function(gene_list, fasta_file, outfile="stdout.fa"){

          ##Read fasta file
          input_sequence = Biostrings::readBStringSet(filepath=fasta_file,use.names = TRUE)
          names(input_sequence) <- base::gsub('[[:space:]| :].*', '', names(input_sequence))


          ID = as.matrix(gene_list)
          #message("no. of ID's: ", nrow(ID))

          ## Find overlap between ID's and fasta file

          overlap <- input_sequence[ID[,1]]

          ## write fasta file

          Biostrings::writeXStringSet(x =overlap, filepath = outfile, format = "fasta")

}


#' Get size of fasta sequence
#' @description Get size (in bp) for each sequence in multi-fasta file.
#'
#' @param fasta_file Either a path or a connection to multi-fasta file.
#'   The input sequence file should have extention .fa or .fasta
#'
#' @return Tab delimited .size file
#' @export
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings width
#' @importFrom utils write.table
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#' @importFrom kableExtra row_spec
#' @importFrom kableExtra column_spec
#' @examples
#' \dontrun{
#'
#'  fa_size(fasta_file="Sc_orf_coding_R64-2-1.fasta")
#'
#' }
fa_size <- function(fasta_file){

          #---- read input genome fasta sequence
          input_sequence <- Biostrings::readBStringSet(filepath=fasta_file,use.names = TRUE)

          #---- extract names and chromosome length
          names(input_sequence) <- base::gsub('[[:space:]| :].*', '', names(input_sequence))

          genome_size_mat <- base::cbind(base::names(input_sequence), Biostrings::width(input_sequence))

          #---- derive output name from input .fa or .fasta file
          output_name <- base::gsub(".fa.*","", basename(fasta_file))

          #---- write output
          utils::write.table(genome_size_mat, file = paste(output_name,".size", sep=""), quote = FALSE,col.names = FALSE, row.names = FALSE,sep = "\t")

          colnames(genome_size_mat) <- c("Seq_id", "Length")

          #message( paste(output_name,".size", sep=""), " file is saved in working directory!")

          temp <- genome_size_mat %>% utils::head(.)

          knitr::kable(temp,"html", align = "l") %>%
                    kableExtra::kable_styling(bootstrap_options = c("striped", "condensed", "responsive"), full_width = F,font_size =14,  stripe_color = "aquamarine3") %>%
                    kableExtra::row_spec(0,bold = TRUE, italic = FALSE, color = "black") %>%
                    kableExtra::column_spec(1:2, bold=FALSE, color="blue")


}

#' Fasta file summary
#' @description Get summary of input multi-fasta file like mean length, median
#'   length, GC content etc
#' @param fasta_file Either a path or a connection to multi-fasta file.
#'   The input sequence file should have extention .fa or .fasta
#'
#' @return A summary table
#' @export
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings width
#' @importFrom gridExtra ttheme_minimal
#' @importFrom gridExtra grid.table
#' @importFrom Biostrings letterFrequency
#' @importFrom tidyr as_tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr summarise
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#' @importFrom dplyr n
#' @importFrom kableExtra row_spec
#' @importFrom kableExtra column_spec
#'
#' @examples
#' \dontrun{
#'
#'  fa_summary(fasta_file="Sc_orf_coding_R64-2-1.fasta")
#' }
fa_summary <- function(fasta_file){

                    # read input sequences as biostring object
                    input_sequence <- Biostrings::readBStringSet(filepath=fasta_file,use.names = TRUE)
                    base::names(input_sequence) <- base::gsub('[[:space:]| :].*', '', base::names(input_sequence))

                    # compute input sequence lengths
                    total_length <-  Biostrings::letterFrequency(input_sequence,letters = "ATGC", collapse=TRUE)
                    gc_content <-  Biostrings::letterFrequency(input_sequence, "GC", collapse=TRUE)


                    summary_length <- cbind(base::names(input_sequence), Biostrings::width(input_sequence)) %>%
                                        tidyr::as_tibble() %>%
                                        dplyr::mutate(V2 =as.numeric(V2)) %>%
                                        dplyr::summarise( num_of_seq=dplyr::n(),min=min(V2), max=max(V2),mean=round(mean(V2),2),
                                                   median= round(median(V2),2)) %>%
                                        dplyr::mutate(percent_gc = round(100*(gc_content/total_length),2))


                    knitr::kable(t(summary_length), col.names = c("Summary"),"html", align = "l") %>%
                              kableExtra::kable_styling(bootstrap_options = c("striped", "condensed", "responsive"), full_width = F,font_size =14,  stripe_color = "aquamarine3") %>%
                              kableExtra::row_spec(0,bold = TRUE, italic = FALSE, color = "black") %>%
                              kableExtra::column_spec(1:2, bold=FALSE, color="blue")



}


#' Get GC percent for multi-fasta file
#' @description for given multi-fasta file get GC percent for each sequence
#' @param fasta_file Either a path or a connection to multi-fasta file. The
#'   input sequence file should have extention .fa or .fasta
#'
#' @return A output file of GC percent and a barplot of GC percent distribution
#'  \code{if the fasta file contains sequences less than 20}
#' @export
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings letterFrequency
#' @importFrom tidyr as_tibble
#' @importFrom dplyr mutate
#' @import ggplot2
#' @importFrom readr write_delim
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#' @importFrom kableExtra row_spec
#' @importFrom kableExtra column_spec
#' @examples
#' \dontrun{
#'
#'  fa_percent_GC(fasta_file="Sc_orf_coding_R64-2-1.fasta")
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
                              dplyr::mutate(percent_gc = round(100*(as.numeric(gc)/as.numeric(sequence_len)),2))


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

                    readr::write_delim(gc_each_seq, path = paste(output_name,"_percentGC.tab"), delim = "\t", col_names = TRUE)

          }

          else{

                    readr::write_delim(gc_each_seq, path = paste(output_name,"_percentGC.tab"), delim = "\t", col_names = TRUE)
          }

          gc_each_seq %>% dplyr::select(c("names","percent_gc")) %>% utils::head() %>%
                    knitr::kable(.,"html", align = "l") %>%
                    kableExtra::kable_styling(bootstrap_options = c("striped", "condensed", "responsive"), full_width = F,font_size =14,  stripe_color = "aquamarine3") %>%
                    kableExtra::row_spec(0,bold = TRUE, italic = FALSE, color = "black") %>%
                    kableExtra::column_spec(1:2, bold=FALSE, color="blue")
}
