#'@title Returns a sequence of a data set of type FASTA
#'@name readSeq
#'
#'@description Given a set of strings in a FASTA file and
#'given an integer X, the sequence X of the set will be returned
#'of data
#'
#'@param arqFasta FASTA type file containing the sequences
#'@param x Integer that indicates which sequence will be returned
#'
#'@details
#'
#'@return Returns the specified string

#'
#' @author Eric Augusto Ito
#'
#' @seealso 
#'
#' @examples
#'
#' 
#' @import seqinr 


lerSeq <- function(arqFasta, x){
	sequencias<-read.fasta(file=arqFasta,forceDNAtolower=FALSE)
	sequencia<-getSequence(sequencias[[x]])
	return(sequencia)
}