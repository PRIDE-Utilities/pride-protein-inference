package uk.ac.ebi.pride.utilities.pia.intermediate;

import uk.ac.ebi.pride.utilities.data.core.DBSequence;

/**
 * An intermediate class, which represents a protein in the intermediate
 * structure.
 * 
 * @author julian
 *
 */
public interface IntermediateProtein {
	
	/**
	 * Getter for an ID
	 * 
	 * @return
	 */
	public Comparable getID();
	
	
	/**
	 * getter for the Protein accession
	 * @return
	 */
	public String getAccession();
	
	
	/**
	 * getter for the protein sequence.
	 * @return null, if no sequence was imported
	 */
	public String getProteinSequence();
	
	
	/**
	 * Getter for the DBSequence representation of this protein.
	 * 
	 * @return
	 */
	public DBSequence getDBSequence();
	
	
	/**
	 * sets the group of this protein
	 * @param group
	 */
	public void setGroup(IntermediateGroup group);
	
	
	/**
	 * returns the group of this protein
	 */
	public IntermediateGroup getGroup();
}
