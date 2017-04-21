#' UbiLength - Ubiquitin interactors of different linear ubiquitin lengths (UbIA-MS dataset)
#'
#' This dataset contains label free quantification (LFQ) data (MaxQuant output \url{http://www.maxquant.org}).
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
#'   \item{Only.identified.by.site}{The protein is only identified by a modification site if marked ('+')}
#'   \item{Reverse}{The protein is identified in the decoy database if marked ('+')}
#'   \item{Potential.contaminant}{The protein is a known contaminant if marked ('+')}
#' }
#' @return A data.frame with the UbiLength dataset
#' @source Zhang, Smits, van Tilburg, et al (2017). An interaction landscape of ubiquitin signaling. Molecular Cell:
#'   \url{http://www.cell.com/molecular-cell/fulltext/S1097-2765(17)30004-7}
#'
"UbiLength"

#' Experimental design of the UbiLength dataset
#'
#' This dataset annotates 12 different samples in 4 conditions and 3 replicates.
#'
#' @format A data.frame with 12 observations and 3 variables:
#' \describe{
#'   \item{label}{Label names}
#'   \item{condition}{Experimental conditions}
#'   \item{replicate}{Replicate number}
#' }
#' @return A data.frame with the UbiLength Experimental design
#' @source Zhang, Smits, van Tilburg, et al (2017). An interaction landscape of ubiquitin signaling. Molecular Cell:
#'   \url{http://www.cell.com/molecular-cell/fulltext/S1097-2765(17)30004-7}
#'
"UbiLength_ExpDesign"

#' GFPip - EED interactors in ESCs (GFP IP - MS dataset)
#'
#' This dataset contains label free quantification (LFQ) and intensity-based absolute quantification (iBAQ) data (MaxQuant output \url{http://www.maxquant.org}).
#'
#' @format A data.frame with 1914 observations and 30 variables:
#' \describe{
#'   \item{Protein.IDs}{Uniprot IDs}
#'   \item{Majority.protein.IDs}{Uniprot IDs of major protein(s) in the protein group}
#'   \item{Protein.names}{Full protein names}
#'   \item{Gene.names}{Gene name}
#'   \item{Fasta.headers}{Header as present in the Uniprot fasta file}
#'   \item{Peptides}{Number of peptides identified for this protein group}
#'   \item{Razor...unique.peptides}{Number of peptides used for the quantification of this protein group}
#'   \item{Unique.peptides}{Number of peptides identified which are unique for this protein group}
#'   \item{Intensity columns (6)}{Raw mass spectrometry intensity, A.U.}
#'   \item{iBAQ columns (6)}{iBAQ normalized mass spectrometry intensity, A.U.}
#'   \item{LFQ.intensity columns (6)}{LFQ normalized mass spectrometry intensity, A.U.}
#'   \item{Only.identified.by.site}{The protein is only identified by a modification site if marked ('+')}
#'   \item{Reverse}{The protein is identified in the decoy database if marked ('+')}
#'   \item{Contaminant}{The protein is a known contaminant if marked ('+')}
#'   \item{id}{The protein group ID}
#' }
#' @return A data.frame with the GFPip dataset
#' @source Kloet et al (2016). The dynamic interactome and genomic targets of Polycomb complexes during stem-cell differentiation. Nature Structural & Molecular Biology:
#'   \url{http://www.nature.com/nsmb/journal/v23/n7/full/nsmb.3248.html}
#'
"GFPip"

#' GFPip_pep - Peptides table of EED interactors in ESCs (GFP IP - MS dataset)
#'
#' This dataset contains the identified peptides (MaxQuant output \url{http://www.maxquant.org}).
#'
#' @format A data.frame with 12750 observations and 29 variables:
#' \describe{
#'   \item{Sequence}{Peptide sequence}
#'   \item{Mass}{Mass of peptide (Da)}
#'   \item{Proteins}{Uniprot IDs of proteins to which the peptide is matched}
#'   \item{Leading.razor.protein}{Uniprot ID of the main protein}
#'   \item{Gene.names}{Gene names of proteins to which the peptide is matched}
#'   \item{Unique..Groups.}{Whether the peptide is unique to the protein group to which the peptide is matched}
#'   \item{Unique..Proteins.}{Whether the peptide is unique to the protein in the complete database}
#'   \item{Charges}{Peptide charge}
#'   \item{PEP}{Posterior Error Probability of the identification}
#'   \item{Score}{Andromeda score}
#'   \item{Reverse}{The peptide is identified in the decoy database if marked ('+')}
#'   \item{Contaminant}{The peptide is a known contaminant if marked ('+')}
#'   \item{id}{The peptide ID}
#'   \item{Protein.group.IDs}{The protein group ID of the protein to which the peptide is matched}
#' }
#' @return A data.frame with the GFPip peptides table
#' @source Kloet et al (2016). The dynamic interactome and genomic targets of Polycomb complexes during stem-cell differentiation. Nature Structural & Molecular Biology:
#'   \url{http://www.nature.com/nsmb/journal/v23/n7/full/nsmb.3248.html}
#'
"GFPip_pep"

#' Experimental design of the GFPip dataset
#'
#' This dataset annotates 6 different samples in 2 conditions and 3 replicates.
#'
#' @format A data.frame with 6 observations and 3 variables:
#' \describe{
#'   \item{label}{Label names}
#'   \item{condition}{Experimental conditions}
#'   \item{replicate}{Replicate number}
#' }
#' @return A data.frame with the GFPip experimental design
#' @source Kloet et al (2016). The dynamic interactome and genomic targets of Polycomb complexes during stem-cell differentiation. Nature Structural & Molecular Biology:
#'   \url{http://www.nature.com/nsmb/journal/v23/n7/full/nsmb.3248.html}
#'
"GFPip_ExpDesign"

#' DiUbi - Ubiquitin interactors for different diubiquitin-linkages (UbIA-MS dataset)
#'
#' This dataset contains label free quantification (LFQ) and intensity-based absolute quantification (iBAQ) data (MaxQuant output \url{http://www.maxquant.org}).
#'
#' @format A data.frame with 4071 observations and 102 variables:
#' \describe{
#'   \item{Protein.IDs}{Uniprot IDs}
#'   \item{Majority.protein.IDs}{Uniprot IDs of major protein(s) in the protein group}
#'   \item{Protein.names}{Full protein names}
#'   \item{Gene.names}{Gene name}
#'   \item{Fasta.headers}{Header as present in the Uniprot fasta file}
#'   \item{Peptides}{Number of peptides identified for this protein group}
#'   \item{Razor...unique.peptides}{Number of peptides used for the quantification of this protein group}
#'   \item{Unique.peptides}{Number of peptides identified which are unique for this protein group}
#'   \item{Intensity columns (30)}{Raw mass spectrometry intensity, A.U.}
#'   \item{iBAQ columns (30)}{iBAQ normalized mass spectrometry intensity, A.U.}
#'   \item{LFQ.intensity columns (30)}{LFQ normalized mass spectrometry intensity, A.U.}
#'   \item{Only.identified.by.site}{The protein is only identified by a modification site if marked ('+')}
#'   \item{Reverse}{The protein is identified in the decoy database if marked ('+')}
#'   \item{Potential.contaminant}{The protein is a known contaminant if marked ('+')}
#'   \item{id}{The protein group ID}
#' }
#' @return A data.frame with the DiUbi dataset
#' @return A data.frame with the GFP Experimental design
#' @source Zhang, Smits, van Tilburg, et al (2017). An interaction landscape of ubiquitin signaling. Molecular Cell:
#'   \url{http://www.cell.com/molecular-cell/fulltext/S1097-2765(17)30004-7}
#'
"DiUbi"

#' Experimental design of the DiUbi dataset
#'
#' This dataset annotates 30 different samples in 10 conditions and 3 replicates.
#'
#' @format A data.frame with 30 observations and 3 variables:
#' \describe{
#'   \item{label}{Label names}
#'   \item{condition}{Experimental conditions}
#'   \item{replicate}{Replicate number}
#' }
#' @return A data.frame with the DiUbi experimental design
#' @source Zhang, Smits, van Tilburg, et al (2017). An interaction landscape of ubiquitin signaling. Molecular Cell:
#'   \url{http://www.cell.com/molecular-cell/fulltext/S1097-2765(17)30004-7}
#'
"DiUbi_ExpDesign"
