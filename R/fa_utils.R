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
#'  faSomeRecords(gene_list=myGenelist, fasta_file="Sc_orf_trans_all_R64-2-1.fasta", outfile="sc_myGenelist.fa")
#'
#' }
faSomeRecords <- function(gene_list, fasta_file, outfile="stdout.fa"){

          ##Read fasta file
          res = Biostrings::readDNAStringSet(filepath=fasta_file,use.names = TRUE)
          names(res) <- base::gsub('[ \\s | :].*', '', names(res))


          ID = as.matrix(gene_list)
          message("no. of ID's: ", nrow(ID))

          ## Find overlap between ID's and fasta file

          overlap <- res[ID[,1]]

          # print(overlap)

          message(cat(outfile)," contains ",length(names(overlap))," Sequences")

          ## write fasta file

          Biostrings::writeXStringSet(x =overlap, filepath = outfile, format = "fasta", use.names=TRUE )

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
#' @examples
#' \dontrun{
#'
#'  faSize(fasta_file="Sc_orf_coding_R64-2-1.fasta")
#'
#' }
faSize <- function(fasta_file){

          #---- read input genome fasta sequence
          input_sequence <- Biostrings::readDNAStringSet(filepath=fasta_file,use.names = TRUE)

          #---- extract names and chromosome length
          names(input_sequence) <- base::gsub(pattern = " .*","",base::names(input_sequence))

          genome_size_mat <- base::cbind(base::names(input_sequence), Biostrings::width(input_sequence))

          #---- derive output name from input .fa or .fasta file
          output_name <- base::gsub(".fa.*","", basename(fasta_file))

          #---- write output
          utils::write.table(genome_size_mat, file = paste(output_name,".size", sep=""), quote = FALSE,col.names = FALSE, row.names = FALSE,sep = "\t")

          colnames(genome_size_mat) <- c("Seq_id", "Length")
          print(knitr::kable(genome_size_mat))

          message( paste(output_name,".size", sep=""), " file is saved in working directory!")

}
