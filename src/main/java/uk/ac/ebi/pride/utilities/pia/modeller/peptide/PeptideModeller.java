package uk.ac.ebi.pride.utilities.pia.modeller.peptide;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateStructure;



/**
 * This class handles the peptides and their inference in a given
 * {@link IntermediateStructure}.
 * 
 * @author julian
 *
 */
public class PeptideModeller {
	
	/** logger for this class */
	private static final Logger logger =  LoggerFactory.getLogger(PeptideModeller.class);
	
	/** the associated intermediate structure */
	private IntermediateStructure intermediateStructure;
	
	
	public PeptideModeller(IntermediateStructure intermediateStructure) {
		this.intermediateStructure = intermediateStructure;
	}
}