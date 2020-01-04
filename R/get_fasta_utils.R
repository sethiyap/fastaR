
#' get fasta sequences for bed file
#'
#' @description Extracts sequences from a fasta file for each of the intervals defined in a BED file
#'
#' @param bedFile A character vector containing path to the file to the bed
#'   file, file format should be standard UCSC bed format with column:
#'   \code{ 1. chromosome-id, 2. start, 3. end, 4. name, 5. score, 6. strand}
#' @param fasta_file Either a path or a connection to reference multi-fasta
#'   file, from which subset of sequences for given input feature is to be
#'   retrieved.
#'   In the sequence header: only string before first space and/or first colon (:)
#'   will be considered for futher processes.
#'   **Important consideration when header have big names.
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
#' @param feature_file A character vector containing path to the file to either
#'   a bed file, file format should be standard UCSC bed format with column:
#'   \code{ 1. chromosome-id, 2. start, 3. end, 4. name, 5. score, 6. strand} OR
#'   gene feature file with extention \code{.gff or .gff3 } ** Chromosome names
#'   should be same as fasta file
#' @param fasta_file Either a path or a connection to reference multi-fasta
#'   file, from which subset of sequences for given input list is to be
#'   retrieved. In the sequence header: only string before first space and/or
#'   first colon (:) will be considered for futher processes. **Important
#'   consideration when header have big names.
#' @param outfile A character vector, determining outfile name, default:
#'   promoter_out
#' @param upstream_bp numeric, base pairs uprstream of start coordinate
#' @param downstream_bp numeric, base pairs downstream of start coordinate
#'
#' @return A output fasta file of promoter sequences and bed file of the
#'   promoter region
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


#' get flanking regions from feature file
#'
#' @description Extracts sequences of one of the flank type described.
#' @param feature_file A character vector containing path to the file to either
#'   a bed file, file format should be standard UCSC bed format with column:
#'   \code{ 1. chromosome-id, 2. start, 3. end, 4. name, 5. score, 6. strand} OR
#'   gene feature file with extention \code{.gff or .gff3 } ** Chromosome names
#'   should be same as fasta file
#' @param width Numeric, width to determine the flank length
#' @param flank_type Numeric,choose region whose sequence (of width length) is
#'   to be fetched. \itemize{ \item 1: sequence upstream of start coordinate
#'   \item 2: sequence downstream of start coordinate \item 3: sequence downstream of end coordinate \item 4: upstream and
#'   downstream of start coordinate \item 5: upstream and downstream of end
#'   coordinate \item 6: upstream and downstream of feature/gene coordinates
#'   \item 7: middle region, ie. width length from start and end coordinate
#'   Start----->ATGCGGATGCGGTC<------End } \code{default: 1}
#' @param outfile A character vector, determining outfile name, \code{default:
#'   flank_out}
#' @param fasta_file Either a path or a connection to reference multi-fasta
#'   file, from which subset of sequences for given input feature is to be
#'   retrieved. In the sequence header: only string before first space and/or
#'   first colon (:) will be considered for futher processes. **Important
#'   consideration when header have big names.
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomicFeatures genes
#' @importFrom stringr str_detect
#' @importFrom rtracklayer import.bed
#' @importFrom GenomicRanges mcols
#' @importFrom Biostrings readDNAStringSet
#' @importFrom IRanges width
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom BSgenome getSeq
#' @importFrom rtracklayer export.bed
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom Biostrings writeXStringSet
#' @importFrom Rsamtools FaFile
#' @importFrom Rsamtools indexFa
#' @importFrom GenomicRanges flank
#' @importFrom IRanges subsetByOverlaps
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#'
#' feature_file <- "Sc_orf_coding_R64-2-1.bed"
#' fasta_file <- "S288C_reference_sequence_R64-2-1_20150113.fa"
#' get_flank_from_feature(feature_file = feature_file, fasta_file = fasta_file, flank_type = 2)
#'
#' }
get_flank_from_feature <- function(feature_file,fasta_file,width=10,flank_type=1,outfile="flank_out"){

          #----- Load Reference feature file

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


          #--- Upstream of start codon
          if(flank_type==1){

                    flank_region <- GenomicRanges::flank(bedfile, start = T,both = F, width=width)

          }

          #--- Downstream of start codon

          if(flank_type==2){

                    flank_region <- GenomicRanges::promoters(bedfile, downstream = width, upstream = 0)
          }

          #--- Downstream of stop codon
          if(flank_type==3){

                    flank_region <- GenomicRanges::flank(bedfile, start = F,both = F, width=width)

          }

          #--- Upstream and Downstream of start codon

          if(flank_type==4){

                    flank_region <- GenomicRanges::flank(bedfile, start = T,both = T, width=width)

          }

          #--- Upstream and Downstream of stop codon
          if(flank_type==5){

                    flank_region <- GenomicRanges::flank(bedfile, start = F,both = T, width=width)

          }

          #--- Upstream and Downstream of CDS region
          if(flank_type==6){

                    flank_region <- bedfile+width

          }

          #--- substract from start and end of CDS region, i.e. how much region within orf from start and end
          #    Start----->---------<------Stop
          if(flank_type==7){

                    flank_region <- bedfile-width

          }

          #---- Get the sequence
          Rsamtools::indexFa(fasta_file)
          genome <- Rsamtools::FaFile(fasta_file)
          dna <- Biostrings::readDNAStringSet(fasta_file)
          names(dna) <- gsub(' .*', '',names(dna))

          #---- Check whether the region boundaries are within genome, remove if out of range

          dd = data.frame(cbind(names(dna), width(dna))) %>%
                    dplyr::mutate(Start=rep(1,length(names(dna)))) %>%
                    dplyr::select(c("X1","Start","X2"))

          colnames(dd)=c("Chr","Start","End")

          dd$Chr <- gsub(' .*', '',dd$Chr)

          dd = GenomicRanges::makeGRangesFromDataFrame(dd)

          flank_region_within_bound <- IRanges::subsetByOverlaps(flank_region,dd,type = "within")

          #---- Get Sequence of within range regions
          flank_seq <- BSgenome::getSeq(dna, flank_region_within_bound)
          names(flank_seq) = flank_region_within_bound$name

          Biostrings::writeXStringSet(flank_seq,filepath=paste(outfile, ".fa",sep=""),append=FALSE,format="fasta" )
          rtracklayer::export.bed(object = flank_region,paste(outfile, ".bed",sep=""))

}

#' get_random_sequences_from_fasta
#'
#' @description Generate n-random sequences of pre-determined length from given input multi-fasta file
#'
#' @param fasta_file Either a path or a connection to reference multi-fasta
#'   file, from which random sequences need to be generated. In the sequence header: only string before first space and/or
#'   first colon (:) will be considered for futher processes. **Important
#'   consideration when header have big names.
#' @param numberOfRandomSequences Numeric, number of random sequences to be generated.
#' @param lengthOfRandomSequence Numeric, length of random sequences.
#' @param outfile A character vector, determining outfile name, \code{default:
#'   randomSeq }
#'
#' @return A bedfile of random regions and fasta file containing random sequences.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' fasta_file <- "S288C_reference_sequence_R64-2-1_20150113.fa"
#'
#' get_random_sequences_from_fasta(fasta_file = fasta_file, numberOfRandomSequences = 100, lengthOfRandomSequence = 500, outfile = "Sc_randomSeq")
#'
#' }
#'
get_random_sequences_from_fasta <- function(fasta_file, numberOfRandomSequences=50, lengthOfRandomSequence=100, outfile="randomSeq"){

          genome.fa <- Biostrings::readDNAStringSet(fasta_file)

          chr_size <- width(genome.fa)
          names(chr_size) <- names(genome.fa)

          span <- 4

          #initialise some vectors for storing random coordinates
          my_random_start  <- vector()
          my_random_end    <- vector()
          my_random_chr    <- vector()
          my_random_strand <- vector()

          set.seed(12345)

          #loop through number of regions
          for(i in 1:numberOfRandomSequences){

                    my_random_chr[i] <- sample(x=names(genome.fa),size=1)
                    my_random_strand[i] <- sample(x=c('-','+'),size=1)
                    my_max <- chr_size[[my_random_chr[i]]]-span
                    my_random_start[i] <- runif(n=1, min=1, max=my_max)
                    my_random_end[i] <- my_random_start[i] + lengthOfRandomSequence
          }

          df <- data.frame(chr=my_random_chr,
                           start=round(my_random_start),
                           end=round(my_random_end),
                           strand=my_random_strand,
                           names=paste("random_seq", seq(1,numberOfRandomSequences, by=1), sep="_"))

          bedfile <- GenomicRanges::makeGRangesFromDataFrame(df,keep.extra.columns = TRUE)

          random_seq = BSgenome::getSeq(genome.fa,bedfile)

          names(random_seq) <- bedfile$names

          Biostrings::writeXStringSet(random_seq,filepath=paste(outfile, ".fa",sep=""),append=FALSE,format="fasta")

          rtracklayer::export.bed(object=bedfile,paste(outfile,".bed",sep=""))

}
