package uk.ac.ebi.pride.utilities.pia.intermediate;

import java.util.List;

import uk.ac.ebi.jmzidml.model.mzidml.AbstractParam;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.pride.utilities.pia.modeller.fdr.FDRComputableByDecoys;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.CvScore;



/**
 * Representation of a peptide spectrum match in the intermediate structure.
 * 
 * @author julian
 *
 */
public abstract class IntermediatePeptideSpectrumMatch implements FDRComputableByDecoys {
	
	/** the connected peptide, if the PSM is already in the intermediate structure */
	private IntermediatePeptide peptide;
	
	/** the decoy status, when overriding the original status */
	private Boolean isDecoy;
	
	/** the calculated FDR value */
	private Double fdrValue;
	
	/** the calculated q-value */
	private Double qValue;
	
	/** the calculated FDR Score value */
	private Double fdrScore;
	
	/** whether this PSM is unique for one protein, does not have to be set by implementation */
	private Boolean isUnique;
	
	
	public IntermediatePeptideSpectrumMatch() {
		this.isDecoy = null;
		this.fdrValue = null;
		this.qValue = null;
		this.fdrScore = null;
		this.isUnique = null;
	}
	
	
	/**
	 * Returns an ID of the peptideSpectrumMatch
	 * 
	 * @return
	 */
	public abstract Comparable getID();
	
	
	/**
	 * Returns the ID of the used controller for importing.
	 * @return
	 */
	public abstract Comparable getControllerID();
	
	
	/**
	 * Returns the score value of the score with the given accession.
	 */
	@Override
	public abstract Double getScore(String scoreAccession);
	
	
	/**
	 * Returns the accessions of the available scores of this PSM, including
	 * calculated scores.
	 */
	public final List<String> getScoreAccessions() {
		List<String> scoreAccessions = getBaseScoreAccessions();
		
		if (getFDR() != null) {
			scoreAccessions.add(CvScore.PSI_PSM_LEVEL_LOCAL_FDR.getAccession());
		}
		
		if (getFDRScore() != null) {
			scoreAccessions.add(CvScore.PSI_PSM_LEVEL_FDRSCORE.getAccession());
		}
		
		if (getQValue() != null) {
			scoreAccessions.add(CvScore.PSI_PSM_LEVEL_Q_VALUE.getAccession());
		}
		
		return scoreAccessions;
	}
	
	
	/**
	 * Returns the basic scores of this PSM. These are the scores saved in the
	 * files and no calculated, additional scores like "FDR Score" or "Q-Value".
	 * 
	 * @return
	 */
	public abstract List<String> getBaseScoreAccessions();
	
	
	/**
	 * Getter for the spectrum id (as in mzIdentML's SpectrumIdentificationResult)
	 */
	public abstract String getSpectrumId();
	
	
	/**
	 * Getter for the PSM's measured experimental m/z
	 */
	public abstract Double getExperimentalMassToCharge();
	
	
	/**
	 * Getter for the mass-delta (NOT m/z) between measured and theoretical
	 * mass value (measured - theoretical)
	 */
	public abstract Double getDeltaMass();
	
	
	/**
	 * Getter for the PSMs' retention time.
	 * <p>
	 * null, if not available
	 */
	public abstract Double getRetentionTime();
	
	
	/**
	 * Getter for the charge.
	 */
	public abstract Integer getCharge();
	
	
	/**
	 * Getter for the number of missed cleavages
	 */
	public abstract Integer getMissedCleavages();
	
	
	/**
	 * Getter for the sequence.
	 */
	public abstract String getSequence();
	
	
	/**
	 * Getter for the modifications 
	 */
	public abstract List<Modification> getModifications();
	
	
	/**
	 * Getter for additional CV and user params (like sourceID, spectrum
	 * title...)
	 * 
	 * @return
	 */
	public abstract List<AbstractParam> getParams();
	
	
	/**
	 * Returns whether the PSM is unique for only one protein (given the
	 * database)
	 */
	public Boolean getIsUnique() {
		if ((isUnique == null) && (peptide != null) && (peptide.getGroup() != null)) {
			// calculate the uniqueness
			switch (peptide.getAllProteins().size()) {
			case 0:
				// uniqueness cannot be calculated, no proteins at all!
				break;
				
			case 1:
				// it is unique
				isUnique = true;
				break;
				
			default:
				// otherwise it is not unique
				isUnique = false;
				break;
			}
		}
		
		return isUnique;
	}
	
	
	/**
	 * Sets the isUnique flag, should only be done by importers.
	 * @param isUnique
	 */
	public void setIsUnique(Boolean isUnique) {
		this.isUnique = isUnique;
	}
	
	
	/**
	 * Getter for the isDecoy flag.
	 * <p>
	 * If the decoy was not set by setIsDecoy, the decoy status of the original
	 * PSM is returned. The original spectrumIdentification is a decoy, if it 
	 * is connected to peptideEvidences / proteins which are only decoys.
	 * 
	 * @return
	 */
	public Boolean getIsDecoy() {
		if (isDecoy != null) {
			return isDecoy;
		} else {
			boolean decoy = true;
			if (peptide != null) {
				for (IntermediateProtein protein : peptide.getAllProteins()) {
					if (protein.getIsDecoy() != null) {
						decoy &= protein.getIsDecoy();
					} else {
						// an un-flagged protein indicates a target
						decoy = false;
					}
					
					if (!decoy) {
						// as soon as it is no longer a decoy, return false
						return false;
					}
				}
			}
			return decoy;
		}
	}
	
	
	/**
	 * Setter for the isDecoy flag.
	 * <p>
	 * This overwrites any decoy settings from the original spectrumIdentification.
	 * This is only for legacy usage of older PIA intermediate files, now the
	 * decoy state is encoded in the protein.
	 */
	@Deprecated
	public void setIsDecoy(Boolean isDecoy) {
		this.isDecoy = isDecoy;
	}
	
	
	/**
	 * Setter for the peptide (called while parsing the PIA XML file).
	 * @param pep
	 */
	public void setPeptide(IntermediatePeptide peptide) {
		this.peptide = peptide;
	}
	
	
	/**
	 * Getter for the peptide. If the peptide is not set while parsing the XML
	 * file, null is returned.
	 * 
	 * @return
	 */
	public IntermediatePeptide getPeptide() {
		return peptide;
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
		if (obj == null || !(obj instanceof IntermediatePeptideSpectrumMatch)) return false;
		
		IntermediatePeptideSpectrumMatch psm = (IntermediatePeptideSpectrumMatch)obj;
		if (!getID().equals(psm.getID())) return false;
		return getControllerID().equals(psm.getControllerID());
	}
	
	
	@Override
	public int hashCode() {
		int result = getID().hashCode();
		result = 31 * getControllerID().hashCode();
		
		return result;
	}
}
