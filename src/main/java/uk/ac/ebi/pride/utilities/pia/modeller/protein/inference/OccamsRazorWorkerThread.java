package uk.ac.ebi.pride.utilities.pia.modeller.protein.inference;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateGroup;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptide;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateProtein;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterUtilities;


public class OccamsRazorWorkerThread extends Thread {
	
	/** the ID of this worker thread */
	private int ID;
	
	/** the caller of this thread */
	private OccamsRazorInference parent;
	
	/** the applied inference filters */
	private List<AbstractFilter> filters;
	
	/** whether modifications are considered while inferring the peptides */
	private boolean considerModifications;
	
	
	
	/** logger for this class */
	private static final Logger logger = Logger.getLogger(OccamsRazorWorkerThread.class);
	
	
	public OccamsRazorWorkerThread(int ID,
			OccamsRazorInference parent,
			List<AbstractFilter> filters,
			boolean considerModifications) {
		this.ID = ID;
		this.parent = parent;
		this.filters = filters;
		this.considerModifications = considerModifications;
		
		this.setName("OccamsRazorWorkerThread-" + this.ID);
	}
	
	
	@Override
	public void run() {
		int treeCount = 0;
		Set<IntermediateGroup> cluster = parent.getNextCluster();
		while (cluster != null) {
			processTree(cluster);
			treeCount++;
			cluster = parent.getNextCluster();
		}
		
		logger.info("<inference worker " + ID + " finished after " + treeCount + " clusters>");
	}
	
	/**
	 * Do the protein inference on the cluster
	 * 
	 * @param cluster
	 */
	private void processTree(Set<IntermediateGroup> cluster) {
		// create the filtered report peptides mapping from the groups' IDs
		Map<Integer, Set<IntermediatePeptide>> groupIdToReportPeptides =
				parent.createClustersFilteredPeptidesMap(cluster, considerModifications);
		
		// the map of actually reported proteins
		List<InferenceProteinGroup> proteins =
				new ArrayList<InferenceProteinGroup>(groupIdToReportPeptides.size());
		
		// maps from the proteinGroup IDs to the peptide keys 
		Map<String, Set<Comparable>> proteinIDsToPeptideIDs = new HashMap<String, Set<Comparable>>();
		
		// create for each group, which has at least one peptide and protein, a protein group
		for (IntermediateGroup group : cluster) {
			
			if ((group.getProteins() == null) ||
					(group.getProteins().size() == 0) ||
					!parent.groupHasReportPeptides(group, groupIdToReportPeptides)) {
				// this group has no proteins or no peptides, skip it
				continue;
			}
			
			// collect the peptides
			Set<IntermediatePeptide> interPeptides = new HashSet<IntermediatePeptide>();
			
			if (groupIdToReportPeptides.containsKey(group.getID())) {
				interPeptides.addAll(groupIdToReportPeptides.get(group.getID()));
			}
			for (IntermediateGroup pepGroup : group.getAllPeptideChildren()) {
				interPeptides.addAll(
						groupIdToReportPeptides.get(pepGroup.getID()));
			}
			
			// look for an existing proteinGroup with same peptides
			boolean addedToExistingGroup = false;
			if ((filters != null ) && (filters.size() > 0)) {
				for (InferenceProteinGroup existingGroup : proteins) {
					
					if (existingGroup.getPeptides().equals(interPeptides)) {
						// add the proteins
						for (IntermediateProtein protein : group.getProteins()) {
							existingGroup.addProtein(protein);
						}
						
						addedToExistingGroup = true;
						break;
					}
				}
			}
			
			// no existing group with same PSMs -> create new group
			if (!addedToExistingGroup) {
				String proteinGroupID = parent.createProteinGroupID(group);
				InferenceProteinGroup proteinGroup = new InferenceProteinGroup(proteinGroupID, considerModifications);
				
				// add the proteins
				for (IntermediateProtein protein : group.getProteins()) {
					proteinGroup.addProtein(protein);
				}
				
				// add the peptides
				proteinGroup.addPeptides(interPeptides);
				
				// and add a peptide mapping
				Set<Comparable> peptides =
						new HashSet<Comparable>(interPeptides.size());
				
				for (IntermediatePeptide peptide : interPeptides) {
					peptides.add(parent.getPeptideKey(peptide, considerModifications));
				}
				
				// calculate the protein's score
				parent.proteinScoring.calculateProteinScore(proteinGroup);
				
				proteinIDsToPeptideIDs.put(proteinGroupID, peptides);
				proteins.add(proteinGroup);
			}
		}
		
		if (proteins.size() < 1) {
			// no proteins could be created (e.g. due to filters?) 
			return;
		}
		
		// remove proteins, not passing the filters
		if ((filters != null ) && (filters.size() > 0)) {
			Iterator<InferenceProteinGroup> proteinIterator = proteins.iterator();
			
			while (proteinIterator.hasNext()) {
				InferenceProteinGroup proteinGroup = proteinIterator.next();
				if (!FilterUtilities.satisfiesFilterList(proteinGroup, filters)) {
					proteinIterator.remove();
				}
			}
		}
		
		// this will be the returned list of proteins
		List<InferenceProteinGroup> reportProteins =
				new ArrayList<InferenceProteinGroup>(proteins.size());
		
		// the still unreported proteins
		HashMap<String, InferenceProteinGroup> unreportedProteins =
				new HashMap<String, InferenceProteinGroup>(proteins.size());
		
		for (InferenceProteinGroup protein : proteins) {
			// add all proteins, which passed the filter, into the map of unreported proteins
			unreportedProteins.put(protein.getID(), protein);
		}
		
		
		// check proteins for sub-proteins and intersections. this cannot be
		// done before, because all proteins have to be built beforehand
		Map<String, Set<InferenceProteinGroup>> proteinIDsToSubproteins =
				new HashMap<String, Set<InferenceProteinGroup>>(groupIdToReportPeptides.size());
		
		Set<InferenceProteinGroup> isSubProtein = new HashSet<InferenceProteinGroup>();
		Set<Comparable> reportedPeptides = new HashSet<Comparable>();
		
		
		
		for (InferenceProteinGroup protein : proteins) {
			Set<Comparable> peptideIDs = proteinIDsToPeptideIDs.get(protein.getID());
			
			Set<InferenceProteinGroup> subProteins = new HashSet<InferenceProteinGroup>();
			proteinIDsToSubproteins.put(protein.getID(), subProteins);
			
			Set<InferenceProteinGroup> intersectingProteins = new HashSet<InferenceProteinGroup>();
			
			boolean thisIsASubProtein = false;
			
			// compare to other proteins
			for (InferenceProteinGroup comparisonProtein : proteins) {
				if (protein == comparisonProtein) {
					// don't compare with same protein
					continue;
				}
				
				Set<Comparable> intersection = new HashSet<Comparable>(
						proteinIDsToPeptideIDs.get(comparisonProtein.getID()));
				intersection.retainAll(peptideIDs);
				
				// there are only intersections now, as the check for samesets was already done
				if (intersection.size() > 0) {
					if (intersection.size() ==
							proteinIDsToPeptideIDs.get(comparisonProtein.getID()).size()) {
						// the complete comparisonProtein is in protein
						subProteins.add(comparisonProtein);
					} else if (intersection.size() == peptideIDs.size()) {
						// the complete proteinID is in subProtID
						isSubProtein.add(protein);
						thisIsASubProtein = true;
					} else if (intersection.size() != peptideIDs.size()) {
						// subProtID intersects proteinID somehow
						intersectingProteins.add(comparisonProtein);
					}
				}
			}
			
			if ((intersectingProteins.size() == 0) &&
					!thisIsASubProtein) {
				// this protein is no subProtein and has no intersections (but
				// maybe subProteins) -> report this protein immediately
				// proteins with intersecting other proteins are reported later
				reportProteins.add(protein);
				reportedPeptides.addAll(proteinIDsToPeptideIDs.get(protein.getID()));
				
				unreportedProteins.remove(protein.getID());
				
				// add the subproteins
				for (InferenceProteinGroup subProtein : subProteins) {
					protein.addSubgroup(subProtein);
					unreportedProteins.remove(subProtein.getID());
				}
			}
		}
		
		// report all the proteins ordered by which explains the most new peptides
		while (unreportedProteins.size() > 0) {
			Set<InferenceProteinGroup> groupsWithMostPeptides = null;
			Set<Comparable> mostCanReport = null;
			int nrMostPeps = -1;
			
			// find the protein group which explains the most not-yet-reported peptides
			for (InferenceProteinGroup protein : unreportedProteins.values()) {
				if (isSubProtein.contains(protein)) {
					// subproteins are reported indirectly with their "parents", skip here
					continue;
				}
				
				Set<Comparable> canReport = proteinIDsToPeptideIDs.get(protein.getID());
				canReport.removeAll(reportedPeptides);
				
				if (canReport.size() > nrMostPeps) {
					groupsWithMostPeptides = new HashSet<InferenceProteinGroup>();
					groupsWithMostPeptides.add(protein);
					nrMostPeps = canReport.size();
					mostCanReport = canReport;
				} else if ((canReport.size() == nrMostPeps) &&
						canReport.equals(mostCanReport)) {
					// another group explains exactly the same peptides -> both are to be reported
					groupsWithMostPeptides.add(protein);
				}
			}
			//Todo: Check probable bug
			for (InferenceProteinGroup protein : groupsWithMostPeptides) {
				reportProteins.add(protein);
				reportedPeptides.addAll(proteinIDsToPeptideIDs.get(protein.getID()));
				
				unreportedProteins.remove(protein.getID());
				
				// add the subproteins
				for (InferenceProteinGroup subProtein : proteinIDsToSubproteins.get(protein.getID())) {
					protein.addSubgroup(subProtein);
					unreportedProteins.remove(subProtein.getID());
				}
			}
		}
		
		if (reportProteins.size() > 0) {
			parent.addToReports(reportProteins);
		}
	}
}
