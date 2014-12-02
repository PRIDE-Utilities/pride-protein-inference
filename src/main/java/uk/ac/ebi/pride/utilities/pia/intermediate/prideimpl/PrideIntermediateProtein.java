package uk.ac.ebi.pride.utilities.pia.intermediate.prideimpl;

import uk.ac.ebi.pride.utilities.data.controller.DataAccessController;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateProtein;


/**
 * An intermediate class, which represents a protein
 * 
 * TODO: implement correct workflow for PRIDE!
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
}
