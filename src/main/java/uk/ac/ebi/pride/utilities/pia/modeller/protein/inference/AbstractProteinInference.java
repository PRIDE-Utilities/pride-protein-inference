package uk.ac.ebi.pride.utilities.pia.modeller.protein.inference;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateGroup;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptide;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateStructure;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterUtilities;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide.PeptideScoring;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.protein.ProteinScoring;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;


/**
 * This abstract class define the main methods for each ProteinInference methods.
 * It works with a controller that do not contain protein groups information. If the
 * data access controller contains ProteinGroup Information the class return the original
 * information from the file. The Protein Inference algorithm also use a set of Filters at
 * at Protein and Peptide level.
 *
 * An {@link AbstractProteinInference} calculates for the PIA groups,
 * which are to be reported and in which way.
 * 
 * @author julian, yperez
 *
 */
public abstract class AbstractProteinInference {
	
	/** the logger for this class */
	private static final Logger logger=  LoggerFactory.getLogger(AbstractProteinInference.class);
	
	/** list of the settings. */
	protected List<AbstractFilter> filters;
	
	/** the currently set peptide scoring */
	protected PeptideScoring peptideScoring;
	
	/** the currently set protein scoring */
	protected ProteinScoring proteinScoring;
	
	/** the number of allowed threads (smaller 1 = all available)*/
	protected int allowedThreads;

	/** the PIA intermediate structure */
	protected IntermediateStructure intermediateStructure;
    
	/** round to two decimals */
	private static final DecimalFormat roundTwoDecimals = new DecimalFormat("#.##");
    
	
    /**
     * The constructor works using a DataAccessController. Any of the current
     * implementations of this class should use this controller to retrieve the
     * psm, peptide and protein information.
     * 
     * @param controller
     * @param filters any filters, which should be applied
     * @param filterPSMsOnImport if true, the PSMs hich fo not pass teh filtering are not imported into the intermediate structure
     * @param nrThreads
     */
	public AbstractProteinInference(IntermediateStructure intermediateStructure,
			PeptideScoring peptideScoring, ProteinScoring proteinScoring,
			List<AbstractFilter> filters, int nrThreads) {
		this.intermediateStructure = intermediateStructure;
		this.filters = (filters == null) ? new ArrayList<AbstractFilter>() : filters;
		this.peptideScoring = peptideScoring;
		this.proteinScoring = proteinScoring;
		this.allowedThreads = nrThreads;
		logger.debug("starting inference with following filters: " + filters);
	}
	
	
    /**
     * calculateInference is the method of the abstract class which allows the
     * class to compute the protein groups and create the List of different
     * groups.
     * 
     * @param considerModifications
     * @param psmSetSettings
     * @return
     */
    public abstract List<InferenceProteinGroup> calculateInference(
    		boolean considerModifications);
    
    
	/**
	 * This method creates a Map from the groups' IDs to the associated
	 * {@link IntermediatePeptide}s, which satisfy the currently set filters.
	 *
	 * @param considerModifications
	 * @return
	 */
	public Map<Integer, Set<IntermediatePeptide>> createFilteredPeptidesMap(boolean considerModifications) {
		Map<Integer, Set<IntermediatePeptide>> groupIdToPeptides =
				new HashMap<Integer, Set<IntermediatePeptide>>(intermediateStructure.getNrGroups() / 2);
		
		for (Set<IntermediateGroup> cluster : intermediateStructure.getClusters().values()) {
			
			Map<Integer, Set<IntermediatePeptide>> clustersGroupIdToPeptides =
					createClustersFilteredPeptidesMap(cluster, considerModifications);
			
			if (clustersGroupIdToPeptides != null) {
				groupIdToPeptides.putAll(clustersGroupIdToPeptides);
			}
		}
		
		return groupIdToPeptides;
	}
	
	
	/**
	 * This method creates a Map from the groups' IDs to the associated
	 * {@link IntermediatePeptide}s of the given cluster, which satisfy the
	 * currently set filters.
	 * 
	 * @param considerModifications
	 * @return
	 */
	public Map<Integer, Set<IntermediatePeptide>> createClustersFilteredPeptidesMap(
			Set<IntermediateGroup> cluster, boolean considerModifications) {
		Map<Integer, Set<IntermediatePeptide>> groupIdToPeptides =
				new HashMap<Integer, Set<IntermediatePeptide>>(intermediateStructure.getNrGroups() / 2);
		
		for (IntermediateGroup group : cluster) {
			
			if ((group.getPeptides() == null) ||
					(group.getPeptides().size() < 1)) {
				// no peptides in the group -> go on
				continue;
			}
			
			Map<Comparable, IntermediatePeptide> groupsPepsMap =
					new HashMap<Comparable, IntermediatePeptide>();
			
			for (IntermediatePeptide pep : group.getPeptides()) {
				if (!considerModifications) {
					// use the same IntermediatePeptide as in the intermediate structure
					groupsPepsMap.put(pep.getSequence(), pep);
					pep.filterPSMs(filters);
				} else {
					for (IntermediatePeptideSpectrumMatch psm : pep.getAllPeptideSpectrumMatches()) {
						if (FilterUtilities.satisfiesFilterList(psm, filters)) {
							// all filters on PSM level are satisfied -> use this PSM
							Comparable pepID = getPSMKey(psm, considerModifications);
							
							// get the peptide of this PSM
							IntermediatePeptide psmsPeptide = groupsPepsMap.get(pepID);
							if (psmsPeptide == null) {
								// no peptide for the pepID in the map yet
								psmsPeptide = new IntermediatePeptide(pep.getSequence());
								groupsPepsMap.put(pepID, psmsPeptide);
							}
							psmsPeptide.addPeptideSpectrumMatch(psm);
						}
					}
				}
			}
			
			Set<IntermediatePeptide> groupsPeptides = new HashSet<IntermediatePeptide>();
			for (IntermediatePeptide pep : groupsPepsMap.values()) {
				if (peptideScoring != null) {
					peptideScoring.calculatePeptideScore(pep);
				}
				
				if ((pep.getPeptideSpectrumMatches().size() > 0) &&
						FilterUtilities.satisfiesFilterList(pep, filters)) {
					// this peptide has PSMs and does satisfy the filters
					groupsPeptides.add(pep);
				}
			}
			groupIdToPeptides.put(group.getID(), groupsPeptides);
		}
		
		return groupIdToPeptides;
	}
	
	
	/**
	 * Tests for the given group, if it has any direct {@link IntermediatePeptide}s,
	 * in the given Map. This Map should by created by 
	 * {@link AbstractProteinInference#createFilteredPeptidesMap(boolean)} prior
	 * to calling this function.
	 */
	public boolean groupHasDirectReportPeptides(IntermediateGroup group,
			Map<Integer, Set<IntermediatePeptide>> groupIdToReportPeptides) {
		Set<IntermediatePeptide> pepSet = groupIdToReportPeptides.get(group.getID());
		return ((pepSet != null) && (pepSet.size() > 0));
	}
	
	
	/**
	 * Tests for the given group, if it has any {@link IntermediatePeptide}s,
	 * whether direct or in the peptideChildren, in the given Map. This Map
	 * should be created by {@link AbstractProteinInference#createFilteredPeptidesMap(boolean)}
	 * prior to calling this function.
	 */
	public boolean groupHasReportPeptides(IntermediateGroup group,
			Map<Integer, Set<IntermediatePeptide>> groupIdToReportPeptides) {
		// check for direct ReportPeptides
		if (groupHasDirectReportPeptides(group, groupIdToReportPeptides)) {
			return true;
		}
		
		// check for ReportPeptides of the children
		for (IntermediateGroup childGroup : group.getAllPeptideChildren()) {
			if (groupHasDirectReportPeptides(childGroup, groupIdToReportPeptides)) {
				return true;
			}
		}
		
		return false;
	}
	
	
	/**
	 * Return a peptide key, depending on whether the modifications are taken
	 * into account or not. The same setting for considerModifications must have
	 * been used while creating the peptide.
	 * 
	 * @param peptide
	 * @param considerModifications
	 * @return
	 */
	static public final Comparable getPeptideKey(IntermediatePeptide peptide, boolean considerModifications) {
		if (considerModifications) {
			// return one of the PSMs' keys
			return getPSMKey(peptide.getAllPeptideSpectrumMatches().get(0), considerModifications);
		} else {
			// just return the sequence
			return peptide.getSequence();
		}
	}
	
	
	/**
	 * Returns a psm key, depending on whether the modifications are taken
	 * into account or not.
	 * 
	 * @param peptide
	 * @param considerModifications
	 * @return
	 */
	static public final Comparable getPSMKey(IntermediatePeptideSpectrumMatch psm, boolean considerModifications) {
		if (considerModifications) {
			StringBuilder sbKey = new StringBuilder(psm.getSequence());
			
			for (Modification mod : psm.getModifications()) {
				sbKey.append('_');
				sbKey.append(mod.getLocation());
				sbKey.append(':');
				sbKey.append(roundTwoDecimals.format(mod.getMonoisotopicMassDelta()));
			}
			
			
			logger.info("Hallo psm " + sbKey.toString());
			return sbKey.toString();
		} else {
			// just return the sequence
			return psm.getSequence();
		}
	}
	
	
	/**
	 * Ceates a unique ID for a ProteinAmbiguityGroup from the group ID
	 * @param group
	 * @return
	 */
	public String createProteinAmbiguityGroupID(IntermediateGroup group) {
		return "PAG_" + group;
	}
	
	
    /**
	 * If polling of inference is performed, return the current state of the
	 * progress (between 0 and 100 progress in percent, &lt;0 inference aborted,
	 * &gt;100 inference is done).
	 * @return
	 */
	public abstract Long getProgressValue();
	
	
	/**
	 * Get the human readable name of the inference method.
	 * @return
	 */
	public abstract String getName();
	
	
	/**
	 * Get the unique machine readable short name of the inference method.
	 * @return
	 */
	public abstract String getShortName();
	
	
	/**
	 * Returns a List of SelectItems representing the available filters for this
	 * inference.
	 *
	 * @return
	 */
	//public abstract List<LabelValueContainer<String>> getFilterTypes();
	
	
	/**
	 * adds a new filter to the inference filters.
	 * @param newFilter
	 * @return
	 */
	public boolean addFilter(AbstractFilter newFilter) {
	    return filters.add(newFilter);
	}
	
	
	/**
	 * Returns a {@link List} of all filter settings for this inference filter.
	 * @return
	 */
	public List<AbstractFilter> getFilters() {
	    return filters;
	}
	
	
	/**
	 * Removes the filter at the given index
	 */
	public AbstractFilter removeFilter(int index) {
		if ((index >= 0) &&
				(index < filters.size())) {
			return filters.remove(index);
		}
		
		return null;
	}
}
