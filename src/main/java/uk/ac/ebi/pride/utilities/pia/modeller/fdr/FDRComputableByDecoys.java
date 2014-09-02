package uk.ac.ebi.pride.utilities.pia.modeller.fdr;


/**
 * Interface used for calculation of FDR estimation by usage of decoys. 
 * 
 * @author julian
 *
 */
public interface FDRComputableByDecoys {
	
	/**
	 * Returns the score value of the score with the given accession.
	 * 
	 * @param scoreAccession
	 * @return
	 */
	public Double getScore(String scoreAccession);
	
	
	/**
	 * Sets the local FDR value.
	 * @return
	 */
	public void setFDR(Double fdr);
	
	
	/**
	 * Gets the local FDR value.
	 * @return
	 */
	public Double getFDR();
	
	
	/**
	 * Getter for the qValue.
	 * @return
	 */
	public Double getQValue();
	
	
	/**
	 * Setter for the qValue.
	 * @return
	 */
	public void setQValue(Double value);
	
	
	/**
	 * Sets the smoothed FDRScore value.
	 * @return
	 */
	public void setFDRScore(Double fdrScore);
	
	
	/**
	 * Gets the smoothed FDRScore value.
	 * @return
	 */
	public Double getFDRScore();
	
	
	/**
	 * Returns true, if the item is a decoy, or false, if not.
	 * @return
	 */
	public Boolean getIsDecoy();
}
