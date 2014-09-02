package uk.ac.ebi.pride.utilities.pia.intermediate;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;


/**
 * This class is a container for the intermediate structure created before
 * PIAs protein inference.
 * 
 * @author julian
 *
 */
public class IntermediateStructure {
	
	/** the clustered groups, the sets of groups are disjoint regarding
	 * peptide-protein connections */
	private Map<Integer, Set<IntermediateGroup>> clusters;
	
	/** the peptides */
	private Set<IntermediatePeptide> peptides;
	
	/** the proteins */
	private Set<IntermediateProtein> proteins;
	
	/** the total number of groups in the structure */
	private Integer nrGroups;
	
	
	public IntermediateStructure() {
		clusters = new HashMap<Integer, Set<IntermediateGroup>>();
		peptides = new HashSet<IntermediatePeptide>();
		proteins = new HashSet<IntermediateProtein>();
		nrGroups = 0;
	}
	
	
	/**
	 * Adds the given cluster to the clusters in the structure. While doing so,
	 * the IDs and treeIDs of the groups are adjusted.
	 * 
	 * @param cluster
	 */
	public void addCluster(Collection<IntermediateGroup> cluster) {
		Integer newTreeID = clusters.size() + 1;
		Set<IntermediateGroup> newGroup = new HashSet<IntermediateGroup>(cluster.size());
		
		for (IntermediateGroup group : cluster) {
			// increase the groups' IDs and adjust the treeIDs
			group.setID(nrGroups++);
			group.setTreeID(newTreeID);
			
			if (group.getPeptides() != null) {
				for (IntermediatePeptide peptide : group.getPeptides()) {
					peptides.add(peptide);
				}
			}
			
			if (group.getProteins() != null) {
				for (IntermediateProtein protein : group.getProteins()) {
					proteins.add(protein);
				}
			}
			
			newGroup.add(group);
		}
		
		clusters.put(newTreeID, newGroup);
	}
	
	
	/**
	 * Returns the number of trees
	 * 
	 * @return
	 */
	public int getNrClusters() {
		return clusters.size();
	}
	
	
	/**
	 * Returns the cluster given by the clusterNr
	 * 
	 * @param clusterNr
	 * @return
	 */
	public Set<IntermediateGroup> getCluster(Integer clusterNr) {
		return clusters.get(clusterNr);
	}
	
	
	/**
	 * Returns the clusters
	 * 
	 * @param clusterNr
	 * @return
	 */
	public Map<Integer, Set<IntermediateGroup>> getClusters() {
		return clusters;
	}
	
	
	/**
	 * Returns the number of total groups
	 * 
	 * @return
	 */
	public int getNrGroups() {
		return nrGroups;
	}
	
	
	/**
	 * Returns the total number of peptides
	 * 
	 * @return
	 */
	public int getNrPeptides() {
		return peptides.size();
	}
	
	
	/**
	 * Returns all intermediate PSMs of the structure
	 * @return
	 */
	public List<IntermediatePeptideSpectrumMatch> getAllIntermediatePSMs() {
		List<IntermediatePeptideSpectrumMatch> psms =
				new ArrayList<IntermediatePeptideSpectrumMatch>(getNrPeptides());
		
		for (IntermediatePeptide pep : peptides) {
			psms.addAll(pep.getAllPeptideSpectrumMatches());
		}
		
		return psms;
	}
}
