#' @title Clean biological secuences
#' @description Curates biological sequences of two restriction enzyme primers or cloning vectors.This cleaning is required for techniques such as cDNA-AFLP.
#' @param SEQs DNAString containing biological sequences that are to be cleaned.
#' @param PrimerF dnastring containing the foward primer/vector sequences to be removed.
#' @param PrimerR dnastring containing the reverse primer/vector sequences to be removed.
#' @return clean biological sequences and visualization of the alignments
#' @author Florencia I Pozzi, Silvina A. Felitti
#' @examples
#' SEQs = DNAString(paste("ACTTTCTGCTGCTTGTGGTCGCAATCAGAGTCCTGATGATGAGTCCTGA",
#' "CCGAACCCTTTTTCTCCGTCATCCGTTGGTCCATGGTACGCAATCAGAG", sep = ""))
#' PrimerF = DNAString("GATGAGTCCTGACCGAA")
#' PrimerR = DNAString("GACTGCGTACCATGC")
#' TwoPrimerRemove (SEQs,PrimerF,PrimerR)
#' @export

TwoPrimerRemove = function (SEQs,PrimerF,PrimerR){
  lengthF = length(PrimerF)
  LENGTHF= (lengthF* 80)/100
  lfF=lengthF-LENGTHF
  lffF=as.integer(lfF)
  m1a = matchPattern(PrimerF, SEQs, max.mismatch=lffF)
  PrimerRRC = reverseComplement(PrimerR)
  nmatchPos1= end(m1a)
  lengthRRC = length(PrimerRRC)
  LENGTHRRC= (lengthRRC* 80)/100
  lfRRC=lengthRRC-LENGTHRRC
  lffRRC=as.integer(lfRRC)
  m2a = matchPattern(PrimerRRC, SEQs, max.mismatch=lffRRC)
  nmatchPos2=start(m2a)
  A=nmatchPos1+1
  B=nmatchPos2-1
  Subseq= DNAStringSet(subseq(SEQs, start=A, end=B))
  print(Subseq)
  fname=tempfile()
  writeXStringSet(Subseq,fname)
  localAlign1 = pairwiseAlignment(PrimerRRC, SEQs, type = "global-local")
  print(localAlign1)
  localAlign2 = pairwiseAlignment(PrimerF, SEQs, type = "global-local")
  print(localAlign2)
}

