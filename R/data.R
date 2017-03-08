#' UbIA_MS - Ubiquitin interactors of different linear ubiquitin lengths
#'
#' This dataset contains label free quantification (LFQ) data (MaxQuant output \url{http://www.maxquant.org})
#'
#' @format A data.frame with 3006 observations and 35 variables:
#' \describe{
#'   \item{Protein.IDs}{Uniprot IDs}
#'   \item{Majority.protein.IDs}{Uniprot IDs of major protein(s) in the protein group}
#'   \item{Protein.names}{Full protein names}
#'   \item{Gene.names}{Gene name}
#'   \item{Fasta.headers}{Header as present in the Uniprot fasta file}
#'   \item{Peptides}{Number of peptides identified for this protein group}
#'   \item{Razor...unique.peptides}{Number of peptides used for the quantification of this protein group}
#'   \item{Unique.peptides}{Number of peptides identified which are unique for this protein group}
#'   \item{Intensity columns (12)}{Raw mass spectrometry intensity, A.U.}
#'   \item{LFQ.intensity columns (12)}{LFQ normalized mass spectrometry intensity, A.U.}
#'   \item{Only.identified.by.site}{The protein is only identified by a modification site if marked ("+")}
#'   \item{Reverse}{The protein is identified in the decoy database if marked ("+")}
#'   \item{Potential.contaminant}{The protein is a known contaminant if marked ("+")}
#' }
#' @source Zhang, Smits, van Tilburg, et al (2017). An interaction landscape of ubiquitin signaling. Molecular Cell
#'   \url{http://www.cell.com/molecular-cell/fulltext/S1097-2765(17)30004-7}
#'
"UbIA_MS"

#' Experimental design of the UbIA_MS dataset
#'
#' This dataset annotates the 12 different labels in 4 samples and 3 replicates.
#'
#' @format A data.frame with 12 observations and 3 variables:
#' \describe{
#'   \item{label}{Uniprot IDs}
#'   \item{sample}{Uniprot IDs of major protein(s) in the protein group}
#'   \item{replicate}{Full protein names}
#' }
#' @source Zhang, Smits, van Tilburg, et al (2017). An interaction landscape of ubiquitin signaling. Molecular Cell
#'   \url{http://www.cell.com/molecular-cell/fulltext/S1097-2765(17)30004-7}
#'
"ExpDesign_UbIA_MS"
