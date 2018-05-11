#' @title Clean biological secuences of one primer
#' @description Cleans biological sequences of one primer.
#' @param SEQ dnastring or file with fasta format containing biological sequences that are to be cleaned.
#' @param PrimerOPR dnastring or file with fasta format containing the reverse primer/vector sequences to be removed.
#' @return clean biological sequences
#' @author Florencia I Pozzi, Gisela Y. Green, Ivana G. Barbona, Gustavo R. Rodriguez, Silvina A. Felitti
#' @examples
#' SEQ4 = DNAString("AATCG")
#' SEQ5 = DNAString("CG")
#' OnePrimerRemove (SEQ4,SEQ5)
#' @import Biostrings
#' @export

OnePrimerRemove = function (SEQ,PrimerOPR){
  PrimerRRA = reverseComplement(PrimerOPR)
  Pos3 = regexpr(PrimerRRA,SEQ)
  nmatchPos3= as.numeric(Pos3)
  Pos4=nchar(SEQ)
  Pos5=Pos4-(Pos4-1)
  nmatchPos6= as.numeric(Pos5)
  A=nmatchPos6
  B=nmatchPos3-1
  Subseq= subseq(SEQ, start=A, end=B)
  return(Subseq)
}

