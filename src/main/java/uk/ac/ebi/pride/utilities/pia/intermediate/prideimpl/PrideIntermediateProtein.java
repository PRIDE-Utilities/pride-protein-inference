package uk.ac.ebi.pride.utilities.pia.intermediate.prideimpl;

import uk.ac.ebi.pride.utilities.data.controller.DataAccessController;
import uk.ac.ebi.pride.utilities.data.core.DBSequence;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateGroup;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateProtein;


/**
 * An intermediate class, which represents a protein
 * 
 * TODO: implement correct workflow for PRIDE!
 * 
 * @author julian
 *
 */
public class PrideIntermediateProtein implements IntermediateProtein {
	
	private DataAccessController controller;
	
	private Comparable proteinID;
	
	private String accession;
	
	/** the connected group of this peptides */
	private IntermediateGroup group;
	
	
	/**
	 * Basic constructor, only initializes the representative
	 *
	 */
	public PrideIntermediateProtein(DataAccessController controller,
			Comparable proteinID) {
		this.controller = controller;
		this.proteinID = proteinID;
		this.accession = controller.getProteinAccession(proteinID);
	}
	
	
	@Override
	public Comparable getID() {
		return accession;
	}
	
	
	@Override
	public String getAccession() {
		return accession;
	}
	
	
	@Override
	public String getProteinSequence() {
		return controller.getProteinSequence(proteinID).getSequence();
	}
	
	
	@Override
	public DBSequence getDBSequence() {
		return controller.getProteinSequence(proteinID);
	}
	
	
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
	public boolean equals(Object obj) {
		if (this == obj) return true;
        if (obj == null || !(obj instanceof IntermediateProtein)) return false;
		
        IntermediateProtein protein = (IntermediateProtein)obj;

        return accession.equals(protein.getAccession()) && !((group != null) ? !group.getID().equals(protein.getGroup().getID()) : (protein.getGroup() != null));
    }
	
	
	@Override
	public int hashCode() {
		int result = accession.hashCode();
        result = 31 * result + ((group != null) ? group.getID().hashCode() : 0);
        return result;
	}
}
