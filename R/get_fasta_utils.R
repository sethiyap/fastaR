
get_fasta_from_bed <- function(bedFile,fasta_file ){

          bedfile = rtracklayer::import.bed(bedFile)
          print(head(bedfile))

          input_sequences = Biostrings::readDNAStringSet(fasta_file,format = "fasta")

          names(input_sequences) <- base::gsub('[ \\s | :].*', '', names(input_sequences))
          print(input_sequences)

          seq = BSgenome::getSeq(input_sequences,bedfile)

          base::names(seq) = bedfile$name
          print(seq)
          Biostrings::writeXStringSet(input_sequences, outfile)
}
