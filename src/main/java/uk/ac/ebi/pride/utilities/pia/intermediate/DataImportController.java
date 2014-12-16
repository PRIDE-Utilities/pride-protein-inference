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
	public void addSpectrumIdentificationsToStructCreator(List<AbstractFilter> filters);
	
	
	/**
	 * Adds the filtered PSMs to the {@link IntermediateStructureCreator}.
	 * Filtering is ok, if the used inference methods are not interfered by it.
	 * <p>
	 * If the controller supports caching of scores, the score with the given
	 * accession will be cached for fast access, which might require slightly
	 * more memory.
	 * 
	 * @param filters
	 * @param cacheScoreAccession
	 */
	public void addSpectrumIdentificationsToStructCreator(List<AbstractFilter> filters, String cacheScoreAccession);
	
	
	/**
	 * Some controllers should be closed after usage.
	 */
	public void close();
}
