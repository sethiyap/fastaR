
#' get fasta sequences for bed file
#'
#' @description Extracts sequences from a fasta file for each of the intervals defined in a BED file
#'
#' @param bedFile A character vector containing path to the file to the bed
#'   file, file format should be standard UCSC bed format with column:
#'   \code{ 1. chromosome-id, 2. start, 3. end, 4. name, 5. score, 6. strand}
#' @param fasta_file Either a path or a connection to reference multi-fasta
#'   file, from which subset of sequences for given input list is to be
#'   retrieved.
#' @param outfile A character vector, determining outfile name, default: out
#'
#' @return A output fasta file of sequences in given bed region
#' @export
#' @importFrom rtracklayer import.bed
#' @importFrom Biostrings readDNAStringSet
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings writeXStringSet
#'
#' @examples
#' \dontrun{
#'
#' get_fasta_from_bed(bedFile="Sc_orf_coding_R64-2-1.bed",
#' fasta_file="S288C_reference_sequence_R64-2-1_20150113.fa",
#' outfile="Sc_subset")
#'
#' }
get_fasta_from_bed <- function(bedFile,fasta_file,outfile="out"){

          bedfile = rtracklayer::import.bed(bedFile)
          print(head(bedfile))

          input_sequences = Biostrings::readDNAStringSet(fasta_file,format = "fasta")

          names(input_sequences) <- base::gsub('[[:space:]| :].*', '', names(input_sequences))

          seq = BSgenome::getSeq(input_sequences,bedfile)

          base::names(seq) = bedfile$name
          print(seq)
          Biostrings::writeXStringSet(seq, filepath = paste0(outfile,"_",length(seq), "sequences.fasta"))
}


#' get fasta sequences of promoters
#' @description Extracts promoter sequences from input bed/gff file
#'
#' @param feature_file A character vector containing path to the file to either a bed
#'   file, file format should be standard UCSC bed format with column:
#'   \code{ 1. chromosome-id, 2. start, 3. end, 4. name, 5. score, 6. strand} OR
#'   gene feature file with extention \code{.gff or .gff3 }
#'   ** Chromosome names should be same as fasta file
#' @param fasta_file Either a path or a connection to reference multi-fasta
#'   file, from which subset of sequences for given input list is to be
#'   retrieved.
#' @param outfile A character vector, determining outfile name, default: promoter_out
#' @param upstream_bp numeric, base pairs uprstream of start coordinate
#' @param downstream_bp numeric, base pairs downstream of start coordinate
#'
#' @return A output fasta file of promoter sequences and bed file of the promoter region
#' @export
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomicFeatures genes
#' @importFrom stringr str_detect
#' @importFrom rtracklayer import.bed
#' @importFrom magrittr set_colnames
#' @importFrom GenomicRanges mcols
#' @importFrom Biostrings readDNAStringSet
#' @importFrom IRanges width
#' @importFrom tidyr as_tibble
#' @importFrom GenomicFeatures promoters
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom BSgenome getSeq
#' @importFrom rtracklayer export.bed
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom Biostrings writeXStringSet
#'
#' @examples
#' \dontrun{
#'
#' get_promoter_from_feature(feature_file="Sc_orf_coding_R64-2-1.bed",
#' fasta_file="S288C_reference_sequence_R64-2-1_20150113.fa",
#' outfile="Sc_promoter", upstream_bp=1000, downstream_bp=100)
#'
#' }
#'
get_promoter_from_feature <- function(feature_file,fasta_file,outfile="promoter_out", upstream_bp=500, downstream_bp=1){
          library(dplyr)

          if(stringr::str_detect(string = as.character(basename(feature_file)),pattern=".gff")){

                    gff = GenomicFeatures::makeTxDbFromGFF(feature_file)
                    bedfile = GenomicFeatures::genes(gff)
                    base::names(bedfile)=NULL
                    gene_name <- names(GenomicRanges::mcols(bedfile)[1])
          }
          if(stringr::str_detect(string = as.character(basename(feature_file)),pattern=".bed")){
                    bedfile = rtracklayer::import.bed(feature_file)
                    gene_name <- names(GenomicRanges::mcols(bedfile)[1])
          }


          input_sequences = Biostrings::readDNAStringSet(fasta_file,format = "fasta")
          names(input_sequences) <- base::gsub('[[:space:]| :].*', '', names(input_sequences))

          sequence_size <- cbind(names(input_sequences), IRanges::width(input_sequences)) %>%
                              tidyr::as_tibble() %>% magrittr::set_colnames(c("V1", "V2")) %>%
                              dplyr::mutate(V2=as.numeric(V2))



          promoter_bed <- GenomicFeatures::promoters(x = bedfile,upstream = upstream_bp,downstream = downstream_bp, use.names = TRUE)

          promoter_tidy <- tidyr::as_tibble(promoter_bed) %>% dplyr::mutate(start= ifelse(start<0, 1, start)) %>%
                                         dplyr::left_join(., sequence_size, by=c("seqnames"="V1")) %>%
                                         dplyr::filter(!(end > V2)) %>%
                                         dplyr::select(c(seqnames, start, end,strand,gene_name))

          # --- make granges of the filtered promoters
          promoter_tidy_filtered <-  GenomicRanges::makeGRangesFromDataFrame(promoter_tidy, keep.extra.columns = T)

          names(GenomicRanges::mcols(promoter_tidy_filtered)) <- "names"

          seq = BSgenome::getSeq(input_sequences,promoter_tidy_filtered)


          base::names(seq) = promoter_tidy_filtered$names
          print(seq)
          rtracklayer::export.bed(object = promoter_tidy_filtered,con = paste0(outfile,"_",length(promoter_bed), "promoters.bed"))
          Biostrings::writeXStringSet(seq, filepath = paste0(outfile,"_",length(seq), "sequences.fasta"))

}
