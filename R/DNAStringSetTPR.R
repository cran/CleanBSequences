#' @title Curing of biological sequences
#' @description Curates biological sequences of two restriction enzyme primers or cloning vectors.This cleaning is required for techniques such as cDNA-AFLP.This cleaning is required for techniques such as cDNA-AFLP.
#' @param SEQs file with fasta format containing biological sequences that are to be cleaned.
#' @param PrimerR dnastring containing the reverse primer/vector sequences to be removed.
#' @param PrimerF dnastring containing the foward primer/vector sequences to be removed.
#' @return clean biological sequences and visualization of the alignments
#' @import pwalign
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings vmatchPattern
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings subseq
#' @importFrom Biostrings writeXStringSet
#' @author Florencia I Pozzi, Silvina A. Felitti
#' @examples
#' SEQs = readDNAStringSet(system.file("sequences","SeqInputTPR.fasta", package = "CleanBSequences"))
#' PrimerR= DNAString ("GACTGCGTACCATGC")
#' PrimerF = DNAString("GATGAGTCCTGACCGAA")
#' DNAStringSetTPR (SEQs,PrimerF,PrimerR)
#' @export

DNAStringSetTPR = function (SEQs,PrimerF,PrimerR){
  lengthF = length(PrimerF)
  LENGTHF= (lengthF* 80)/100
  lfF=lengthF-LENGTHF
  lffF=as.integer(lfF)
  m1a = vmatchPattern(PrimerF, SEQs, max.mismatch=lffF)
  nmatchPos1= end(m1a)
  PrimerRRC = reverseComplement(PrimerR)
  lengthRRC = length(PrimerRRC)
  LENGTHRRC= (lengthRRC* 80)/100
  lfRRC=lengthRRC-LENGTHRRC
  lffRRC=as.integer(lfRRC)
  m2a = vmatchPattern(PrimerRRC, SEQs, max.mismatch=lffRRC)
  nmatchPos2=start(m2a)
  A=nmatchPos1+1
  B=nmatchPos2-1
  C = as.integer(A)
  D = as.integer(B)
  Subseq= DNAStringSet(subseq(SEQs, start=C, end=D))
  print(Subseq)
  fname=tempfile()
  writeXStringSet(Subseq,fname)
  localAlign = pairwiseAlignment(SEQs, PrimerRRC, type = "local")
  print(localAlign)
  localAlign2 = pairwiseAlignment(SEQs,PrimerF, type = "local")
  print(localAlign2)
}

