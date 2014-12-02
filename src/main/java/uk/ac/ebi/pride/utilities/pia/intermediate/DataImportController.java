package uk.ac.ebi.pride.utilities.pia.intermediate;

import java.util.List;

import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;


/**
 * An interface for the data import controllers for PIA input files.
 * 
 * @author julian
 *
 */
public interface DataImportController {
	
	/**
	 * Returns the ID of the controller.
	 * 
	 * @return
	 */
	public Comparable getID();
	
	
	/**
	 * Returns the file name which is handled by this controller
	 * 
	 * @return
	 */
	public String getInputFileName();
	
	
	/**
	 * Adds the filtered PSMs to the {@link IntermediateStructureCreator}.
	 * Filtering is ok, if the used inference methods are not interfered by it.
	 * 
	 * @param structCreator
	 */
	public void addSpectrumIdentificationsToStructCreator(IntermediateStructureCreator structCreator, List<AbstractFilter> filters);
	
	
	/**
	 * Some controllers should be closed after usage.
	 */
	public void close();
}
