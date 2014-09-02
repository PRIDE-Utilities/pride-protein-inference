package uk.ac.ebi.pride.utilities.pia.intermediate;

import uk.ac.ebi.pride.utilities.data.core.SpectrumIdentification;
import uk.ac.ebi.pride.utilities.pia.modeller.fdr.FDRComputableByDecoys;


/**
 * Representation of a peptide spectrum match in the intermediate structure.
 * 
 * @author julian
 *
 */
public interface IntermediatePeptideSpectrumMatch extends FDRComputableByDecoys {
	
	/**
	 * Returns an ID of the peptideSpectrumMatch
	 * 
	 * @return
	 */
	public Comparable getID();
	
	
	/**
	 * Returns the ID of the used controller for importing.
	 * @return
	 */
	public Comparable getControllerID();
	
	
	/**
	 * Getter for the isDecoy flag.
	 * <p>
	 * If the decoy was not set by setIsDecoy, the decoy status of the original
	 * spectrumIdentification is returned. The original spectrumIdentification
	 * is a decoy, if it contains PeptideEvidences with decoys only.
	 * 
	 * @return
	 */
	@Override
	public Boolean getIsDecoy();
	
	
	/**
	 * Setter for the isDecoy flag.
	 * <p>
	 * This overwrites any decoy settings from the original spectrumIdentification.
	 */
	public void setIsDecoy(Boolean isDecoy);
	
	
	/**
	 * Returns the actual spectrum identification.
	 * 
	 * @return
	 */
	public SpectrumIdentification getSpectrumIdentification();
}
