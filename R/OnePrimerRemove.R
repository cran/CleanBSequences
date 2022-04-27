#' @title Clean biological secuences
#' @description Curates biological sequences of primer reverse.This cleaning is required for techniques such as cDNA-AFLP.
#' @param SEQs dnastring or file with fasta format containing biological sequences that are to be cleaned.
#' @param PrimerR dnastring or file with fasta format containing the reverse primer/vector sequences to be removed.
#' @return clean biological sequences and visualization of the alignments
#' @import Biostrings
#' @author Florencia I Pozzi and Silvina A. Felitti
#' @examples
#' SEQs = DNAString(paste("GCCTCGCCTCCCTCTTTGATCAGCTTCGCATATCAGGCAACAGCTCAATTT",
#' "GGTACTTGTTCAAATAAGCATTTAGACCATCTGTTCCAAGAACCTTTGCAATCTT",
#' "CACAAGGTGGTCATGGTACGCAGTC", sep=""))
#' PrimerR= DNAString("GACTGCGTACCATGC")
#' OnePrimerRemove (SEQs,PrimerR)
#' @export

OnePrimerRemove = function (SEQs,PrimerR){
  PrimerRRC = reverseComplement(PrimerR)
  lengthRRC = length(PrimerRRC)
  LENGTHRRC= (lengthRRC* 80)/100
  lfRRC=lengthRRC-LENGTHRRC
  lffRRC=as.integer(lfRRC)
  ma = matchPattern(PrimerRRC, SEQs, max.mismatch=lffRRC)
  nmatchPos2=start(ma)
  A=1
  B=nmatchPos2-1
  Subseq= DNAStringSet(subseq(SEQs, start=A, end=B))
  print(Subseq)
  fname=tempfile()
  writeXStringSet(Subseq,fname)
  localAlign = pairwiseAlignment(PrimerRRC, SEQs, type = "global-local")
  print(localAlign)
}

