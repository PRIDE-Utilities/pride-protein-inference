package uk.ac.ebi.pride.utilities.pia.modeller.protein.inference;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateGroup;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptide;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateProtein;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateStructure;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterUtilities;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide.PeptideScoring;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.protein.ProteinScoring;


/**
 * This inference reports all the PIA {@link IntermediateGroup}s as one protein.
 * <p>
 * This is similar to distinguish proteins simply by their peptides and report
 * every possible set and subset.
 * 
 * @author julian
 *
 */
public class ReportAllInference extends AbstractProteinInference {
	
	/** the logger for this class */
	private static final Logger logger= Logger.getLogger(ReportAllInference.class);
	
	/** the human readable name of this filter */
	protected static final String name = "Report All";
	
	/** the machine readable name of the filter */
	protected static final String shortName = "inference_report_all";
	
	/** the progress of the inference */
	private Double progress;
	
	
	public ReportAllInference(IntermediateStructure intermediateStructure,
			PeptideScoring peptideScoring, ProteinScoring proteinScoring,
			List<AbstractFilter> filters, int nrThreads) {
		super(intermediateStructure, peptideScoring, proteinScoring, filters, nrThreads);
		this.progress = 0.0;
	}
	
	/*
	@Override
	public List<LabelValueContainer<String>> getFilterTypes() {
		List<LabelValueContainer<String>> filters = new ArrayList<LabelValueContainer<String>>();
		
		filters.add(new LabelValueContainer<String>(null, "--- PSM ---"));
		for (Map.Entry<String, String>  scoreIt
				: getAvailableScoreShorts().entrySet()) {
			String[] filterNames = PSMScoreFilter.getShortAndFilteringName(
					scoreIt.getKey(), scoreIt.getValue());
			
			if (filterNames != null) {
				filters.add(new LabelValueContainer<String>(filterNames[0], filterNames[1]));
			}
		}
		filters.add(new LabelValueContainer<String>(NrPSMsPerPSMSetFilter.shortName(), NrPSMsPerPSMSetFilter.filteringName()));
		
		return filters;
		return null;
	}
	*/
	
	@Override
	public List<InferenceProteinGroup> calculateInference(boolean considerModifications) {
		progress = 0.0;
		logger.info("calculateInference started...");
		/*
		logger.info("scoring: " + getScoring().getName() + " with " + 
				getScoring().getScoreSetting().getValue() + ", " +
				getScoring().getPSMForScoringSetting().getValue());
		*/
		
		// all the PSMs of the groups, including the PSMs in groups' children
		Map<Integer, Set<IntermediatePeptide>> groupsAllPeptides =
				new HashMap<Integer, Set<IntermediatePeptide>>(intermediateStructure.getNrGroups());
		
		// maps from the groups' IDs to the groups  with equal PSMs after filtering
		Map<Integer, Set<IntermediateGroup>> sameSets = null;
		
		// the finally returned list of protein groups
		List<InferenceProteinGroup> proteinGroups = new ArrayList<InferenceProteinGroup>();
		
		for (Set<IntermediateGroup> cluster : intermediateStructure.getClusters().values()) {
			
			// maps from the groups' IDs to the peptides, which should be reported
			Map<Integer, Set<IntermediatePeptide>> groupIdToReportPeptides =
				createClustersFilteredPeptidesMap(cluster, considerModifications);
			
			Double progressStep = 89.0 / intermediateStructure.getNrClusters() / cluster.size();
			Set<IntermediateGroup> clusterReportGroups = new HashSet<IntermediateGroup>(cluster.size());
			
			// put every group with direct accessions into the report map map
			for (IntermediateGroup group : cluster) {
				if (((group.getProteins() != null) && (group.getProteins().size() > 0)) &&
						groupHasReportPeptides(group, groupIdToReportPeptides)) {
					// report this group
					clusterReportGroups.add(group);
					
					// get the peptides of this group
					Set<IntermediatePeptide> allPeptidesSet = new HashSet<IntermediatePeptide>();
					groupsAllPeptides.put(group.getID(), allPeptidesSet);
					
					// add the direct peptides
					if (groupIdToReportPeptides.containsKey(group.getID())) {
						for (IntermediatePeptide peptide : groupIdToReportPeptides.get(group.getID())) {
							allPeptidesSet.add(peptide);
						}
					}
					
					// add childrens' peptides
					for (IntermediateGroup pepGroup : group.getAllPeptideChildren()) {
						if (groupIdToReportPeptides.containsKey(pepGroup.getID())) {
							for (IntermediatePeptide peptide : groupIdToReportPeptides.get(pepGroup.getID())) {
								allPeptidesSet.add(peptide);
							}
						}
					}
				}
				
				progress += progressStep;
			}
			
			// check for sameSets (if there were active filters)
			if ((filters != null ) && (filters.size() > 0)) {
				
				sameSets = new HashMap<Integer, Set<IntermediateGroup>>(groupsAllPeptides.size());
				Set<Integer> newReportGroupIDs = new HashSet<Integer>(clusterReportGroups.size());
				
				for (Map.Entry<Integer, Set<IntermediatePeptide>> gIt : groupsAllPeptides.entrySet()) {
					// every group gets a sameSet
					Set<IntermediateGroup> sameSet = sameSets.get(gIt.getKey()); 
					if (sameSet == null) {
						sameSet = new HashSet<IntermediateGroup>();
						sameSets.put(gIt.getKey(), sameSet);
					}
					
					// check against the other report groups
					for (IntermediateGroup checkGroup : clusterReportGroups) {
						if (gIt.getKey() == checkGroup.getID()) {
							// don't check against self
							continue;
						}
						
						if (gIt.getValue().equals(groupsAllPeptides.get(checkGroup.getID()))) {
							// ReportPeptides are the same in checkSet and grIt
							sameSet.add(checkGroup);
							
							// if checkID's group had a sameSet before, merge the sameSets
							Set<IntermediateGroup> checkSameSet = sameSets.get(checkGroup.getID());
							if (checkSameSet != null) {
								sameSet.addAll(checkSameSet);
							}
							sameSets.put(checkGroup.getID(), sameSet);
						}
					}
					
					// check, if any of the sameSet is already in the newReportGroups 
					boolean anySameInReportGroups = false;
					
					for (IntermediateGroup sameGroup : sameSet) {
						if (newReportGroupIDs.contains(sameGroup.getID())) {
							anySameInReportGroups = true;
							break;
						}
					}
					
					if (!anySameInReportGroups) {
						// no sameGroup in reportGroups yet, put this one in
						newReportGroupIDs.add(gIt.getKey());
					}
				}
				
				Set<IntermediateGroup> newReportGroups = new HashSet<IntermediateGroup>(newReportGroupIDs.size());
				for (IntermediateGroup group : clusterReportGroups) {
					if (newReportGroupIDs.contains(group.getID())) {
						newReportGroups.add(group);
					}
				}
				clusterReportGroups = newReportGroups;
			}
			
			// now create the proteins from the groups, which are in clusterReportGroups
			for (IntermediateGroup group : clusterReportGroups) {
				InferenceProteinGroup proteinGroup =
						new InferenceProteinGroup(createProteinGroupID(group), considerModifications);
				
				// add the proteins to the group
				for (IntermediateProtein interProt : group.getProteins()) {
					proteinGroup.addProtein(interProt);
				}
				
				if (sameSets != null) {
					for (IntermediateGroup sameGroup : sameSets.get(group.getID())) {
						if (sameGroup.getID() == group.getID()) {
							continue;
						}
						
						for (IntermediateProtein interProt : sameGroup.getProteins()) {
							proteinGroup.addProtein(interProt);
						}
					}
				}
				
				for (IntermediatePeptide interPeptide : groupsAllPeptides.get(group.getID())) {
					// peptides are already filtered in createClustersFilteredPeptidesMap
					proteinGroup.addPeptide(interPeptide);
				}
				
				if (proteinScoring != null) {
					proteinScoring.calculateProteinScore(proteinGroup);
				}
				
				if (FilterUtilities.satisfiesFilterList(proteinGroup, filters)) {
					// add only proteinGroups, which satisfy the filtering
					proteinGroups.add(proteinGroup);
				}
			}
			
			progress += 10.0 / intermediateStructure.getNrClusters();
		}
		
		logger.info("calculateInference done.");
		progress = 100.0;
		return proteinGroups;
	}
	
	
	@Override
	public String getName() {
		return name;
	}
	
	
	@Override
	public String getShortName() {
		return shortName;
	}

	@Override
	public Long getProgressValue() {
		return progress.longValue();
	}
}