#' @title Curing of biological sequences
#' @description Curates biological sequences of primer reverse.This cleaning is required for techniques such as cDNA-AFLP.
#' @param SEQs file with fasta format containing biological sequences that are to be cleaned.
#' @param PrimerR dnastring containing the reverse primer/vector sequences to be removed.
#' @return clean biological sequences and visualization of the alignments
#' @import Biostrings
#' @author Florencia I Pozzi, Silvina A. Felitti
#' @examples
#' SEQs = readDNAStringSet(system.file("sequences","SeqInputOPR.fasta", package = "CleanBSequences"))
#' PrimerR= DNAString ("GACTGCGTACCATGC")
#' DNAStringSetOPR (SEQs,PrimerR)
#' @export

DNAStringSetOPR = function (SEQs,PrimerR){
  PrimerRRC = reverseComplement(PrimerR)
  lengthRRC = length(PrimerRRC)
  LENGTHRRC= (lengthRRC* 80)/100
  lfRRC=lengthRRC-LENGTHRRC
  lffRRC=as.integer(lfRRC)
  ma = vmatchPattern(PrimerRRC, SEQs, max.mismatch=lffRRC)
  nmatchPos2=start(ma)
  A=1
  B=nmatchPos2-1
  C = as.integer(B)
  Subseq= DNAStringSet(subseq(SEQs, start=A, end=C))
  print(Subseq)
  fname=tempfile()
  writeXStringSet(Subseq,fname)
  localAlign = pairwiseAlignment(SEQs, PrimerRRC, type = "local")
  print(localAlign)
}
