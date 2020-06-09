
#' get fasta sequences from bed file
#'
#' @description Extracts sequences from a fasta file for each of the intervals defined in a BED file
#'
#' @param bedFile Either a path or a connection to the bed
#'   file, file format should be standard UCSC bed format with column:
#'   \code{ 1. chromosome-id, 2. start, 3. end, 4. name, 5. score, 6. strand}
#'   ** The file should be without column names.
#' @param fasta_file Either a path or a connection to reference multi-fasta
#'   file, from which subset of sequences for given input feature is to be
#'   retrieved.
#'   In the sequence header: only string before first space and/or first colon (:)
#'   will be considered for further processes.
#'   **Important consideration when header have long names.
#' @param outfile Logical, to return output sequences as a output multi-fasta file or not, default: FALSE
#'
#' @return sequences in given bed region
#' @export
#' @importFrom rtracklayer import.bed
#' @importFrom Biostrings readDNAStringSet
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings writeXStringSet
#'
#' @examples
#' \dontrun{
#'
#' bed_file_in <- system.file("exdata","Sc_ref_genes.bed", package = "fastaR")
#' ref_fasta <- system.file("exdata", "Sc_ref_genome.fasta", package = "fastaR")
#' fastaR::get_fasta_from_bed(bedFile=bed_file_in, fasta_file=ref_fasta, write_output=FALSE)
#'
#' }
get_fasta_from_bed <- function(bedFile,fasta_file,write_output=FALSE){

          bedfile = rtracklayer::import.bed(bedFile)
          print(head(bedfile))

          input_sequences = Biostrings::readDNAStringSet(fasta_file,format = "fasta")

          names(input_sequences) <- gsub('[[:space:]| :].*', '', names(input_sequences))

          seq = BSgenome::getSeq(input_sequences,bedfile)

          names(seq) = bedfile$name

          if(write_output==TRUE){

                    output_name <- gsub(".bed.*","", basename(bedFile))
                    Biostrings::writeXStringSet(seq, filepath = paste0(output_name,"_",length(seq), "sequences.fasta"))
          } else{
                  return(seq)
          }

}


#' get fasta sequences of promoters
#' @description Extracts promoter sequences from input bed/gff file
#' @param feature_file Either a path or a connection to either
#'   a bed file, file format should be standard UCSC bed format with column:
#'   \code{ 1. chromosome-id, 2. start, 3. end, 4. name, 5. score, 6. strand} OR
#'   gene feature file with extention \code{.gff or .gff3 } ** Chromosome names
#'   should be same as fasta file
#' @param fasta_file Either a path or a connection to reference multi-fasta
#'   file, from which subset of sequences for given input list is to be
#'   retrieved. In the sequence header: only string before first space and/or
#'   first colon (:) will be considered for further processes. **Important
#'   consideration when header have long names.
#' @param upstream_bp numeric, base pairs upstream of start coordinate
#' @param downstream_bp numeric, base pairs downstream of start coordinate
#' @param write_outputfasta Logical, to return promoter sequences as a output multi-fasta file, Default: FALSE
#' @param write_promoterbed Logical, to return promoter regions as a output bed file, Default: FALSE
#' @param upstream_bp numeric, base pairs window upstream of start coordinate, Default: 500
#' @param downstream_bp numeric, base pairs window downstream of start coordinate, Default: 1
#' @return sequences of promoters and bed file of the promoter region
#' @export
#' @rdname get_promoter_from_feature
#' @importFrom stringr str_detect
#' @importFrom GenomicFeatures makeTxDbFromGFF genes promoters
#' @importFrom GenomicRanges mcols makeGRangesFromDataFrame
#' @importFrom rtracklayer import.bed export.bed
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom dplyr bind_cols mutate left_join filter select all_of
#' @importFrom IRanges width
#' @importFrom tidyr as_tibble
#' @importFrom BSgenome getSeq
#' @examples
#' \dontrun{
#'
#' feature_file_in <- system.file("exdata","Sc_ref_genes.gff", package = "fastaR")
#' ref_fasta <- system.file("exdata", "Sc_ref_genome.fasta", package = "fastaR")
#' fastaR::get_promoter_from_feature(feature_file=feature_file_in, fasta_file=ref_fasta,
#' write_outputfasta= FALSE, write_promoterbed= FALSE, upstream_bp= 1000, downstream_bp= 100)
#'
#' }
get_promoter_from_feature <- function(feature_file,fasta_file,write_outputfasta=FALSE, write_promoterbed=FALSE, upstream_bp=500, downstream_bp=1){


          if(stringr::str_detect(string = as.character(basename(feature_file)),pattern=".gff")){

                    gff = GenomicFeatures::makeTxDbFromGFF(feature_file, format = "gff3")
                    bedfile = GenomicFeatures::genes(gff)
                    names(bedfile)=NULL
                    gene_name <- names(GenomicRanges::mcols(bedfile)[1])
          }
          if(stringr::str_detect(string = as.character(basename(feature_file)),pattern=".bed")){
                    bedfile = rtracklayer::import.bed(feature_file)
                    gene_name <- names(GenomicRanges::mcols(bedfile)[1])
          }


          input_sequences = Biostrings::readDNAStringSet(fasta_file,format = "fasta")
          names(input_sequences) <- gsub('[[:space:]| :].*', '', names(input_sequences))

          sequence_size <- dplyr::bind_cols(V1=names(input_sequences),V2= IRanges::width(input_sequences)) %>%
                              dplyr::mutate(V2=as.numeric(V2))

          promoter_bed <- GenomicFeatures::promoters(x = bedfile,upstream = upstream_bp,downstream = downstream_bp, use.names = TRUE)

          promoter_tidy <- tidyr::as_tibble(promoter_bed) %>% dplyr::mutate(start= ifelse(start<0, 1, start)) %>%
                                         dplyr::left_join(., sequence_size, by=c("seqnames"="V1")) %>%
                                         dplyr::filter(!(end > V2)) %>%
                                         dplyr::select(c(seqnames, start, end,strand,dplyr::all_of(gene_name)))

          # --- make granges of the filtered promoters

          promoter_tidy_filtered <-  GenomicRanges::makeGRangesFromDataFrame(promoter_tidy, keep.extra.columns = T)

          names(GenomicRanges::mcols(promoter_tidy_filtered)) <- "names"

          seq = BSgenome::getSeq(input_sequences,promoter_tidy_filtered)

          names(seq) = promoter_tidy_filtered$names

          output_name <- gsub('[.].*','', basename(feature_file))

          if(write_promoterbed==TRUE){

                    rtracklayer::export.bed(object = promoter_tidy_filtered,con = paste0(output_name,"_",length(promoter_bed), "promoters.bed"))
          } else {
                    print(head(promoter_tidy_filtered))
          }
          if(write_outputfasta==TRUE){

                    Biostrings::writeXStringSet(seq, filepath = paste0(output_name,"_",length(seq), "promoter_sequences.fasta"))
          } else {
                    return(seq)
          }
}


#' get flanking regions from feature file
#'
#' @description Extracts sequences of one of the flank type described.
#' @param feature_file Either a path or a connection to either
#'   a bed file, file format should be standard UCSC bed format with column:
#'   \code{ 1. chromosome-id, 2. start, 3. end, 4. name, 5. score, 6. strand} OR
#'   gene feature file with extention \code{.gff or .gff3 } ** Chromosome names
#'   should be same as fasta file.
#' @param width Numeric, width to determine the flank length, Default: 10
#' @param flank_type Numeric,choose region whose sequence (of width length) is
#'   to be fetched. \itemize{ \item 1: sequence upstream of start coordinate
#'   \item 2: sequence downstream of start coordinate \item 3: sequence downstream of end coordinate \item 4: upstream and
#'   downstream of start coordinate \item 5: upstream and downstream of end
#'   coordinate \item 6: upstream and downstream of feature/gene coordinates
#'   \item 7: middle region, ie. width length from start and end coordinate
#'   Start----->ATGCGGATGCGGTC<------End } Default: 1
#' @param fasta_file Either a path or a connection to reference multi-fasta
#'   file, from which subset of sequences for given input feature is to be
#'   retrieved. In the sequence header: only string before first space and/or
#'   first colon (:) will be considered for further processes. **Important
#'   consideration when header have long names.
#' @param write_flankbed Logical, to return flank region as a output bed file, Default: FALSE
#' @param write_outputfasta Logical, to return flank sequences as a output multi-fasta file, Default: FALSE
#' @param outfile character vector, defining output file name, Default: 'flank_out'
#' @importFrom stringr str_detect
#' @importFrom GenomicFeatures makeTxDbFromGFF genes
#' @importFrom rtracklayer import.bed export.bed
#' @importFrom GenomicRanges flank promoters makeGRangesFromDataFrame
#' @importFrom Rsamtools indexFa FaFile
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom dplyr bind_cols mutate select
#' @importFrom IRanges width subsetByOverlaps
#' @importFrom BSgenome getSeq
#' @return fasta sequence and ranges of the flank region
#' @export
#'
#' @examples
#' \dontrun{
#'
#' feature_file_in <- system.file("exdata","Sc_ref_genes.gff", package = "fastaR")
#' ref_fasta <- system.file("exdata", "Sc_ref_genome.fasta", package = "fastaR")
#' fastaR::get_flank_from_feature(feature_file = feature_file_in, fasta_file = ref_fasta, flank_type = 2)
#'
#' }
get_flank_from_feature <- function(feature_file,fasta_file,width=10,flank_type=1,write_flankbed=FALSE, write_outputfasta=FALSE, outfile="flank_out"){

          #----- Load Reference feature file

          if(stringr::str_detect(string = as.character(basename(feature_file)),pattern=".gff")){

                    gff = GenomicFeatures::makeTxDbFromGFF(feature_file, format = "gff")
                    bedfile = GenomicFeatures::genes(gff)
                    names(bedfile)=NULL
          }
          if(stringr::str_detect(string = as.character(basename(feature_file)),pattern=".bed")){
                    bedfile = rtracklayer::import.bed(feature_file)
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

          dd = dplyr::bind_cols(Chr=names(dna), End=IRanges::width(dna)) %>%
                    dplyr::mutate(Start=rep(1,length(names(dna)))) %>%
                    dplyr::select(c(Chr,Start,End))

          #colnames(dd)=c("Chr","Start","End")

          dd$Chr <- gsub(' .*', '',dd$Chr)

          dd = GenomicRanges::makeGRangesFromDataFrame(dd)

          flank_region_within_bound <- IRanges::subsetByOverlaps(flank_region,dd,type = "within")

          #---- Get Sequence of within range regions
          flank_seq <- BSgenome::getSeq(dna, flank_region_within_bound)

          if(stringr::str_detect(string = as.character(basename(feature_file)),pattern=".bed")){
                    names(flank_seq) <-  flank_region_within_bound$name
          }
          if(stringr::str_detect(string = as.character(basename(feature_file)),pattern=".gff")){
                    names(flank_seq) <-  flank_region_within_bound$gene_id
          }

          # write output or return
          if(write_flankbed==TRUE){
                    rtracklayer::export.bed(object = flank_region,paste(outfile, ".bed",sep=""))
          } else {
                    print(head(flank_region))
          }
          if(write_outputfasta==TRUE){

                    Biostrings::writeXStringSet(flank_seq,filepath=paste(outfile, ".fa",sep=""),append=FALSE,format="fasta" )
          } else {
                    return(flank_seq)
          }
}

#' get_random_sequences_from_fasta
#'
#' @description Generate n-random sequences of pre-determined length from given input multi-fasta file
#'
#' @param fasta_file Either a path or a connection to reference multi-fasta
#'   file, from which subset of sequences for given input feature is to be
#'   retrieved. In the sequence header: only string before first space and/or
#'   first colon (:) will be considered for further processes. **Important
#'   consideration when header have long names.
#' @param numberOfRandomSequences Numeric, number of random sequences to be generated.
#' @param lengthOfRandomSequence Numeric, length of random sequences.
#' @param write_bed Logical, to return flank region as a output bed file, Default: FALSE
#' @param write_outputfasta Logical, to return flank sequences as a output multi-fasta file, Default: FALSE
#' @param outfile character vector, defining output file name, Default: 'flank_out'
#' @return A ranges of random regions and respective random sequences.
#' @rdname get_random_sequences_from_fasta
#' @export
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom IRanges width
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom BSgenome getSeq
#' @importFrom rtracklayer export.bed
#' @examples
#' \dontrun{
#'
#' ref_fasta <- system.file("exdata", "Sc_ref_genome.fasta", package = "fastaR")
#'
#' fastaR::get_random_sequences_from_fasta(fasta_file = ref_fasta,
#'                                 numberOfRandomSequences = 100,
#'                                 lengthOfRandomSequence = 500)
#'
#' }
#'
get_random_sequences_from_fasta <- function(fasta_file, numberOfRandomSequences=50, lengthOfRandomSequence=100, write_bed=FALSE,write_outputfasta=FALSE,outfile="randomSeq"){

          genome.fa <- Biostrings::readDNAStringSet(fasta_file)

          chr_size <- IRanges::width(genome.fa)
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

          # write output or return
          if(write_bed==TRUE){
                    rtracklayer::export.bed(object=bedfile,paste(outfile,".bed",sep=""))
          } else {
                    print(head(bedfile))
          }
          if(write_outputfasta==TRUE){

                    Biostrings::writeXStringSet(random_seq,filepath=paste(outfile, ".fa",sep=""),append=FALSE,format="fasta")
          } else {
                    return(random_seq)
          }
}
