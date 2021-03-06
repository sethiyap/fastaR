% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fa_utils.R
\name{fa_size}
\alias{fa_size}
\title{Get size of fasta sequence}
\usage{
fa_size(fasta_file)
}
\arguments{
\item{fasta_file}{Either a path or a connection to multi-fasta file.
The input sequence file should have extention .fa or .fasta
In the sequence header: only string before first space and/or first colon (:) will be considered for futher processes.
**Important consideration when header have big names.}
}
\value{
A tibble of gene_id and gene length
}
\description{
Get size (in bp) for each sequence in multi-fasta file.
}
\examples{
\dontrun{

 ref_fasta <- system.file("exdata", "Sc_nucl_R64-2-1.fasta", package = "fastaR")

 fastaR::fa_size(fasta_file=ref_fasta)

}
}
