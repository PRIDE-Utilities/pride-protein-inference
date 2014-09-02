package uk.ac.ebi.pride.utilities.pia.intermediate;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;


/**
 * This thread builds up the intermediate structure given the peptide to
 * dbSequence mapping. 
 * 
 * @author julian
 *
 */
public class IntermediateStructureCreatorWorkerThread extends Thread {
	
	/** the ID of this worker thread */
	private int ID;
	
	/** the parent/creator of this thread */
	private IntermediateStructureCreator parent;
	
	
	/** logger for this class */
	private static final Logger logger = Logger.getLogger(IntermediateStructureCreatorWorkerThread.class);
	
	
	public IntermediateStructureCreatorWorkerThread(int ID, IntermediateStructureCreator parent) {
		this.ID = ID;
		this.parent = parent;
		
		this.setName(getClass().getCanonicalName()+"-"+ID);
	}
	
	
	@Override
	public void run() {
		int workedClusters = 0;
		Map<Comparable, Set<Comparable>> cluster;
		
		// get the next available cluster from the parent
		cluster = parent.getNextCluster();
		while (cluster != null) {
			// the created groups of this cluster
			Map<Integer, IntermediateGroup> subGroups =
					new HashMap<Integer, IntermediateGroup>();
			
			// mapping from the dbSequence IDs to the intermediateProteins
			Map<Comparable, IntermediateProtein> dbSeqsToProteins =
					new HashMap<Comparable, IntermediateProtein>();
			
			for (Map.Entry<Comparable, Set<Comparable>> pepIt : cluster.entrySet()) {
				
				IntermediatePeptide pep = parent.getPeptide(pepIt.getKey());
				
				Set<IntermediateProtein> proteins =
						new HashSet<IntermediateProtein>(pepIt.getValue().size());
				for (Comparable proteinID : pepIt.getValue()) {
					proteins.add(parent.getProtein(proteinID));
				}
				
				insertIntoMap(pep, proteins, subGroups);
			}
			
			// put the subGroups as new tree into the intermediateStructure
			parent.addCluster(subGroups.values());
			workedClusters++;
			
			// get the next cluster information
			cluster = parent.getNextCluster();
		}
		
		logger.info("<thread " + ID + " has no more work after " +
				workedClusters + " clusters> ");
	}
	
	
	/**
	 * Inserts the given peptide with its accessions into the intermediate
	 * format, which is then build up by the groups in subGroups.
	 * 
	 * @param peptide
	 * @param dbSequences
	 * @param subGroups
	 */
	public void insertIntoMap(IntermediatePeptide peptide,
			Set<IntermediateProtein> proteins,
			Map<Integer, IntermediateGroup> subGroups) {
		
		Map<Integer, Set<IntermediateProtein>> groupToProteinsMap;	// the proteins, grouped by their groups' IDs
		
		// group the dbSequences by the intermediate groups they are in
		groupToProteinsMap =
				new HashMap<Integer, Set<IntermediateProtein>>(subGroups.size());
		for (IntermediateProtein protein : proteins) {
			Integer id = (protein.getGroup() == null) ? -1 : protein.getGroup().getID();
			
			Set<IntermediateProtein> groupsProteins = groupToProteinsMap.get(id);
			if (groupsProteins == null) {
				groupsProteins = new HashSet<IntermediateProtein>();
				groupToProteinsMap.put(id, groupsProteins);
			}
			
			groupsProteins.add(protein);
		}
		
		if ((groupToProteinsMap.size() == 1) && (groupToProteinsMap.containsKey(-1))) {
			// all dbSequences are not yet assigned to any group
		    //  => assign all to a new group
			
			// create the new group
			IntermediateGroup group = new IntermediateGroup(subGroups.size()+1);
			subGroups.put(group.getID(), group);
			
			// connect peptide and group
			connectPeptideToGroup(peptide, group);
			
			// add all accessions to this new group
			for (IntermediateProtein protein : proteins) {
				connectProteinToGroup(protein, group);
			}
		} else {
			if (groupToProteinsMap.size() == 1) {
				// all accessions have the same group (but are already assigned)
				// get group of the accessions
				
				// there is only one id, so get this group
				IntermediateGroup group = 
						subGroups.get(groupToProteinsMap.keySet().iterator().next());
				
				// look, if the group of the dbSequences has any other dbSequence(s)
				if (groupHasOnlyTheseDirectProteins(group, proteins)) {
					// the group of the accessions has only the accessions to be
					// assigned to the peptide
			        //   => add the peptide to the group
					connectPeptideToGroup(peptide, group);
				} else {
			        // the group of the accessions has NOT only the accessions
					// to be assigned to the peptide
			        //   => create a group with the peptide and the accessions
					//      group and move the accessions there
					
			        // create the new group
					IntermediateGroup betweenGroup =
							new IntermediateGroup(subGroups.size()+1);
					subGroups.put(betweenGroup.getID(), betweenGroup);
					
			        // add group to peptide and vice versa
					connectPeptideToGroup(peptide, betweenGroup);

			        // add the old group to new group as child
					connectGroups(betweenGroup, group);
					
					for (IntermediateProtein protein : groupToProteinsMap.get(group.getID())) {
						connectProteinToGroup(protein, betweenGroup);
					}
				}
			} else {
				// (accGrouped.size() != 1)
				// the accessions are in different groups / not yet assigned
				
				Set<Integer> remainingGroups = new HashSet<Integer>();
				Set<Integer> subTreeSet = getSubtreeGroups(proteins,
						remainingGroups, subGroups);
				
				if ((remainingGroups.size() == 0) &&
						(((subTreeSet.size() == 1) && !subTreeSet.contains(-1L)) ||			// either there is only one group (and it's not -1, the unassigned)
								((subTreeSet.size() == 2) && subTreeSet.contains(-1L)))) {	// or there are 2 groups and one of it are the unassigned accessions
					// the already assigned accessions build up a whole subtree
					
			        // get the group building up the subtree (the one with the assigned)
					IntermediateGroup group = null;
					for (Integer id : subTreeSet) {
						if (id > 0) {
							group = subGroups.get(id);
						}
					}
					
					if (group != null) {
						if (groupToProteinsMap.containsKey(-1L)) {
							// we have some unassigned accessions as well
							
							// create a between group
							IntermediateGroup betweenGroup =
									new IntermediateGroup(subGroups.size()+1);
							subGroups.put(betweenGroup.getID(), betweenGroup);
							
							// add the unassigned accessions to the between group
							for (IntermediateProtein protein : groupToProteinsMap.get(-1)) {
								connectProteinToGroup(protein, betweenGroup);
							}
							
							// add the new group as child to the group
							connectGroups(group, betweenGroup);
							
							group = betweenGroup;
						}
						
						// add group to peptide and vice versa
						connectPeptideToGroup(peptide, group);
					} else {
						logger.error("There should have been a group for the accessions!");
					}
				} else {
					// can't say much about the constellation of groups
					
					// create new group for peptide
					IntermediateGroup pepGroup =
							new IntermediateGroup(subGroups.size()+1);
					subGroups.put(pepGroup.getID(), pepGroup);
					
					// add group to peptide and vice versa
					connectPeptideToGroup(peptide, pepGroup);
					
					// add the new group to all subTreeSet-groups as new child
					
					for (Integer subTreeId : subTreeSet) {
						if (subTreeId == -1) {
							// if we have unassigned accessions, add them
							// directly to the pepGroup
							for (IntermediateProtein protein : groupToProteinsMap.get(-1)) {
								connectProteinToGroup(protein, pepGroup);
							}
						} else {
							// add the pepGroup to the other (not unassigned
							// accessions) groups
							connectGroups(subGroups.get(subTreeId), pepGroup);
						}
					}
					
					// move the accessions of the remainingGroups into a new
					// group, pointing to the pepGroup and the old group
					for (Integer remGroupId : remainingGroups) {
						// the "remaining" group
						IntermediateGroup group = subGroups.get(remGroupId);
						
						// create an between group
						IntermediateGroup betweenGroup =
								new IntermediateGroup(subGroups.size()+1);
						subGroups.put(betweenGroup.getID(), betweenGroup);
						
						// connect the between group to the remaining and the pepGroup
						connectGroups(betweenGroup, group);
						connectGroups(betweenGroup, pepGroup);
						
						// move the accessions to the between group
						for (IntermediateProtein protein : groupToProteinsMap.get(remGroupId)) {
							connectProteinToGroup(protein, betweenGroup);
						}
					}
				}
			}
		}
	}
	
	
	/**
	 * Connects the given peptide with the given group. A warning is given, if
	 * the peptide was connected to another group.
	 *  
	 * @param peptide
	 * @param group
	 */
	private void connectPeptideToGroup(IntermediatePeptide peptide, IntermediateGroup group) {
		if (peptide.getGroup() != null) {
			logger.warn("peptide " + peptide.getSequence() +
					" was already connected to a group!");
		}
		peptide.setGroup(group);
		group.addPeptide(peptide);
	}
	
	
	/**
	 * Connects the two given groups.
	 * 
	 * @param parent
	 * @param child
	 */
	private void connectGroups(IntermediateGroup parent, IntermediateGroup child) {
		parent.addChild(child);
		child.addParent(parent);
	}
	
	/**
	 * Returns false if the group has any dbSequences which are not in the given
	 * set or any parent dbSequences.
	 * 
	 * @param group
	 * @param dbSeqeunces
	 * @return
	 */
	private boolean groupHasOnlyTheseDirectProteins(IntermediateGroup group,
			Set<IntermediateProtein> proteins) {
		
		// if the group has parent-groups, it has (by default) other dbSequences
		if ((group.getParents() != null) && (group.getParents().size() > 0)) {
			return false;
		}
		
		// check if all of the group's proteins are in the given set
        return !((group.getProteins() != null) &&
                !proteins.containsAll(group.getProteins()));

    }
	
	
	/**
	 * Returns false if the group has any dbSequences in allDbSequences which
	 * are not in the given set.
	 * 
	 * @param group
	 * @param proteins
	 * @return
	 */
	private boolean groupHasOnlyTheseInAllProteins(IntermediateGroup group,
			Set<IntermediateProtein> proteins) {
		// iterate through the map of the group's allAccessions and look for
		// accessions not in accessions
		for (IntermediateProtein protein : group.getAllProteins()) {
			if (!proteins.contains(protein)) {
				return false;
			}
		}
		
		return true;
	}
	
	
	/**
	 * Connects the given protein with the given group.
	 * <p>
	 * If the accession was connected to another group before, it will be
	 * removed there.
	 * 
	 * @param accession
	 * @param group
	 */
	private void connectProteinToGroup(IntermediateProtein protein, IntermediateGroup group) {
		if ((protein.getGroup() != null) &&
				!protein.getGroup().equals(group)) {
			protein.getGroup().removeProtein(protein);
		}
		
		protein.setGroup(group);
		group.addProtein(protein);
	}
	
	
	/**
	 * Builds a set with the group IDs which are the lowest
	 * (closest to the peptides) in the tree to build up the subtrees with only
	 * the accessions in the given set.
	 * <p>
	 * The variable remainingSet will be filled with the group's IDs of
	 * dbSequences, which are not satisfied with the built tree, so their groups
	 * have other dbSequences which are not in the set of dbSequences.
	 * 
	 * @param proteins
	 * @param remainingSet
	 * @param subGroups
	 * @return
	 */
	private Set<Integer> getSubtreeGroups(Set<IntermediateProtein> proteins,
			Set<Integer> remainingSet, Map<Integer, IntermediateGroup> subGroups) {
		Set<Integer> subTreeSet = new HashSet<Integer>();
		IntermediateGroup mostGroup = null;
		
		// look for unassigned dbSequences
		for (IntermediateProtein protein : proteins) {
			if (protein.getGroup() == null) {
				// if we have (at least) one unassigned accession, insert the -1 (unassigned) group
				subTreeSet.add(-1);
				break;
			}
		}
		
		do {
			int nrMostGroupDbSequences = 0;
			mostGroup = null;
			
			// get the group (in subGroups) with the most dbSequences in its
			// allAccessions, which are also in the accessions map
			for (IntermediateGroup group : subGroups.values()) {
				if (groupHasOnlyTheseInAllProteins(group, proteins)) {
					Set<IntermediateProtein> allProteins = group.getAllProteins(); 
					if ((allProteins.size() > nrMostGroupDbSequences) ||
							((mostGroup == null) && (allProteins.size() > 0))) {
						nrMostGroupDbSequences = allProteins.size();
						mostGroup = group;
					}
				}
			}
			
			// remove the accessions of the found group from the accessions set
			if (mostGroup != null) {
				for (IntermediateProtein protein : mostGroup.getAllProteins()) {
					proteins.remove(protein);
				}
				
				// and add the found ID to the subTreeSet
				subTreeSet.add(mostGroup.getID());
			}
		} while (mostGroup != null);
		
		remainingSet.clear();
		for (IntermediateProtein protein : proteins) {
			if (protein.getGroup() != null) {
				// if there is still an accession in the map, which has an
				// assigned group, so there were accessions in the group, which
				// were not in the accessions map -> put it into the remainingSet
				remainingSet.add(protein.getGroup().getID());
			}
		}
		
		return subTreeSet;
	}
}
