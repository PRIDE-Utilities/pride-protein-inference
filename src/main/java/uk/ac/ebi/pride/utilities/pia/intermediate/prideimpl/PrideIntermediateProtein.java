package uk.ac.ebi.pride.utilities.pia.intermediate.prideimpl;

import uk.ac.ebi.pride.utilities.data.controller.DataAccessController;
import uk.ac.ebi.pride.utilities.data.core.Protein;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateProtein;


/**
 * An intermediate class, which represents a protein
 * 
 * @author julian
 *
 */
public class PrideIntermediateProtein extends IntermediateProtein {
	
	private DataAccessController controller;
	
	private Comparable proteinID;
	
	private String accession;
	
	
	/**
	 * Basic constructor, only initializes the representative
	 *  
	 * @param sequence
	 */
	public PrideIntermediateProtein(DataAccessController controller,
			Comparable proteinID) {
		super();
		this.controller = controller;
		this.proteinID = proteinID;
		this.accession = controller.getProteinAccession(proteinID);
	}
	
	
	@Override
	public String getAccession() {
		return accession;
	}
	
	
	@Override
	public String getProteinSequence() {
		return controller.getProteinSequence(proteinID).getSequence();
	}
	
	
	/**
	 * Getter for the PRIDE protein ID
	 * @return
	 */
	public Comparable getPrideProteinID() {
		return proteinID;
	}
	
	
	/**
	 * Getter for the complete PRIDE protein
	 * @return
	 */
	public Protein getPrideProtein() {
		return controller.getProteinById(proteinID);
	}
}
