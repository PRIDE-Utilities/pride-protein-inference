package uk.ac.ebi.pride.utilities.pia.intermediate;

import uk.ac.ebi.pride.utilities.pia.modeller.fdr.FDRComputableByDecoys;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.CvScore;



/**
 * An intermediate class, which represents a protein in the intermediate
 * structure.
 * 
 * @author julian
 *
 */
public abstract class IntermediateProtein implements FDRComputableByDecoys {
	
	/** the connected group of this peptides */
	private IntermediateGroup group;
	
	/** the set isDecoy flag */
	private Boolean isDecoy;
	
	/** the calculated FDR value */
	private Double fdrValue;
	
	/** the calculated q-value */
	private Double qValue;
	
	/** the calculated FDR Score value */
	private Double fdrScore;
	
	/** the protein score (if it is calculatet */
	private Double proteinScore;
	
	
	public IntermediateProtein() {
		this.group = null;
		this.isDecoy = null;
		this.fdrValue = null;
		this.qValue = null;
		this.fdrScore = null;
		this.proteinScore = null;
	}
	
	
	/**
	 * getter for the Protein accession
	 * @return
	 */
	public abstract String getAccession();
	
	
	/**
	 * getter for the protein sequence.
	 * @return null, if no sequence was imported
	 */
	public abstract String getProteinSequence();
	
	
	/**
	 * sets the group of this protein
	 * @param group
	 */
	public void setGroup(IntermediateGroup group) {
		this.group = group;
	}
	
	
	/**
	 * returns the group of this protein
	 */
	public IntermediateGroup getGroup() {
		return group;
	}
	
	
	@Override
	public Boolean getIsDecoy() {
		if (isDecoy == null) {
			/* TODO: as soon as yasset has implemented this, we will have it
			for (Peptide pep : controller.getProteinById(proteinID).getPeptides()) {
				// check if peptide is decoy....
				
			}
			*/
			return null;
		} else {
			// the isDecoy flag was overwritten, return the set value
			return isDecoy;
		}
	}
	
	
	/**
	 * Setter for the isDecoy flag, may overwrite the initial setting.
	 * @param isDecoy
	 */
	public void setIsDecoy(Boolean isDecoy) {
		this.isDecoy = isDecoy;
	}
	
	
	@Override
	public void setFDR(Double fdr) {
		this.fdrValue = fdr;
	}
	
	
	@Override
	public Double getFDR() {
		return fdrValue;
	}
	
	
	@Override
	public void setQValue(Double value) {
		this.qValue = value;
	}
	
	
	@Override
	public Double getQValue() {
		return qValue;
	}
	
	
	@Override
	public void setFDRScore(Double fdrScore) {
		this.fdrScore = fdrScore;
	}
	
	
	@Override
	public Double getFDRScore() {
		return fdrScore;
	}
	
	
	@Override
	public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || !(obj instanceof IntermediateProtein)) return false;

        IntermediateProtein protein = (IntermediateProtein) obj;

        return getAccession().equals(protein.getAccession()) && !((group != null) ? !group.getID().equals(protein.getGroup().getID()) : (protein.getGroup() != null));
    }
	
	
	@Override
	public int hashCode() {
		int result = getAccession().hashCode();
        result = 31 * result + ((group != null) ? group.getID().hashCode() : 0);
        return result;
	}
	
	
	/**
	 * Returns the score value of the score with the given accession.
	 * <p>
	 * The protein only has one score, i.e. the protein score, and this is
	 * null, if not yet calculated.
	 * 
	 * @param scoreAccession
	 * @return
	 */
	@Override
	public Double getScore(String scoreAccession) {
		if (CvScore.PSI_PIA_PROTEIN_SCORE.getAccession().equals(scoreAccession)) {
			return proteinScore;
		}
		
		return null;
	}
	
	
	/**
	 * Returns the proteinScore or null, if it is not yet calculated.
	 * @return
	 */
	public Double getProteinScore() {
		return proteinScore;
	}
	
	
	/**
	 * Setter for the protein score
	 * @param score
	 */
	public void setProteinScore(Double score) {
		this.proteinScore = score;
	}
}
