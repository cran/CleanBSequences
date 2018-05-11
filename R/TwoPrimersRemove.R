#' @title Clean biological secuences of two primers
#' @description Clean biological sequences of two restriction enzyme primers or cloning vectors.
#' @param SEQs file with fasta format containing biological sequences that are to be cleaned.
#' @param PrimerF dnastring or file with fasta format containing the foward primer/vector sequences to be removed.
#' @param PrimerR dnastring or file with fasta format containing the reverse primer/vector sequences to be removed.
#' @return clean biological sequences
#' @author Florencia I Pozzi, Gisela Y. Green, Ivana G. Barbona, Gustavo R. Rodriguez, Silvina A. Felitti
#' @examples
#' SEQ1 = DNAString("AATCG")
#' SEQ2 = DNAString("TT")
#' SEQ3 = DNAString("CG")
#' TwoPrimerRemove (SEQ1,SEQ2,SEQ3)
#' @import Biostrings
#' @export

TwoPrimerRemove = function (SEQs,PrimerF,PrimerR){
  PrimerFC = complement(PrimerF)
  PrimerRRC = reverseComplement(PrimerR)
  Pos = regexpr(PrimerFC,SEQs)
  Pos1 = regexpr(PrimerRRC,SEQs)
  length= attr(Pos,"match.length")
  nlength=as.numeric(length)
  nmatchPos=as.numeric(Pos)
  nmatchPos1= as.numeric(Pos1)
  A=nmatchPos+nlength
  B=nmatchPos1-1
  Subseq= subseq(SEQs, start=A, end=B)
  return(Subseq)
}
