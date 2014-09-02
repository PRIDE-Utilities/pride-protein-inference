package uk.ac.ebi.pride.utilities.pia.intermediate;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

/**
 * This class creates the intermediate structure needed for fast data access
 * during protein inference.
 * 
 * @author julian
 *
 */
public class IntermediateStructureCreator {
	
	/** the logger for this class */
	private static final Logger logger= Logger.getLogger(IntermediateStructureCreator.class);
	
	/** mapping from the peptide ID to intermediatePeptide
	 * TODO: we could decide in this class, whether a peptide is defined by the sequence only or also by the mods
	 **/
	private Map<Comparable, IntermediatePeptide> peptides;
	
	/** mapping from protein ID to the protein */
	private Map<Comparable, IntermediateProtein> proteins;
	
	
	
	/** mapping from the protein IDs to connected peptides' IDs **/
	private Map<Comparable, Set<Comparable>> proteinsToPeptidesMapping;
	
	/** mapping from the peptide IDs to connected proteins' IDs */
	private Map<Comparable, Set<Comparable>> peptidesToProteinsMapping;
	
	
	
	/** iterates over the clustered list of peptide -> proteins mapping */
	private ListIterator<Map<Comparable, Set<Comparable>>> clusterIterator;
	
	
	
	
	/** the created intermediate structure */
	private IntermediateStructure intermediateStructure;
	
	/** the maximal number of used threads */
	private int numberThreads;
	
	
	public IntermediateStructureCreator(int threads) {
		this.peptides = new HashMap<Comparable, IntermediatePeptide>();
		this.proteins = new HashMap<Comparable, IntermediateProtein>();
		this.proteinsToPeptidesMapping =
				new HashMap<Comparable, Set<Comparable>>();
		this.peptidesToProteinsMapping = new HashMap<Comparable, Set<Comparable>>();
		
		this.clusterIterator = null;
		this.intermediateStructure = null;
		
		this.numberThreads = threads;
	}
	
	
	/**
	 * returns true, if the peptides map already contains a peptide with the
	 * given ID
	 * @param peptideID
	 * @return
	 */
	public boolean peptidesContains(Comparable peptideID) {
		return peptides.containsKey(peptideID);
	}
	
	
	/**
	 * adds the given peptide to the peptides map
	 * 
	 * @param peptide
	 * @return any previous peptide with the same ID or null
	 */
	public IntermediatePeptide addPeptide(IntermediatePeptide peptide) {
		return peptides.put(peptide.getID(), peptide);
	}
	
	
	/**
	 * getter for the peptide from the map with the given ID
	 * @param pepId
	 * @return
	 */
	public IntermediatePeptide getPeptide(Comparable pepId) {
		return peptides.get(pepId);
	}
	
	
	/**
	 * returns true, if the proteins map already contains a protein with the
	 * given ID
	 * @param proteinID
	 * @return
	 */
	public boolean proteinsContains(Comparable proteinID) {
		return proteins.containsKey(proteinID);
	}
	
	
	/**
	 * adds the given protein to the proteins map
	 * 
	 * @param protein
	 * @return any previous protein with the same ID or null
	 */
	public IntermediateProtein addProtein(IntermediateProtein protein) {
		return proteins.put(protein.getID(), protein);
	}
	
	
	/**
	 * getter for the protein from the map with the given ID
	 * 
	 * @param proteinID
	 * @return 
	 */
	public IntermediateProtein getProtein(Comparable proteinID) {
		return proteins.get(proteinID);
	}
	
	
	/**
	 * Connects the peptide with the protein. Both have to be already in the
	 * peptides respectively proteins map.
	 * 
	 * @param peptideID
	 * @param proteinID
	 * @return
	 */
	public void addPeptideToProteinConnection(Comparable peptideID, Comparable proteinID) {
		Set<Comparable> protIDs = peptidesToProteinsMapping.get(peptideID);
		if (protIDs == null) {
			protIDs = new HashSet<Comparable>();
			peptidesToProteinsMapping.put(peptideID, protIDs);
		}
		protIDs.add(proteinID);
		
		Set<Comparable> pepIDs = proteinsToPeptidesMapping.get(proteinID);
		if (pepIDs == null) {
			pepIDs = new HashSet<Comparable>();
			proteinsToPeptidesMapping.put(proteinID, pepIDs);
		}
		pepIDs.add(peptideID);
	}
	
	
	public int getNrPeptides() {
		return peptides.size();
	}
	
	
	public int getNrProteins() {
		return proteins.size();
	}
	
	
	public int getNrSpectrumIdentifications() {
		int nrSpectrumIdentifications = 0;
		for (IntermediatePeptide pep : peptides.values()) {
			nrSpectrumIdentifications += pep.getAllPeptideSpectrumMatches().size();
		}
		return nrSpectrumIdentifications;
	}
	
	
	/**
	 * After the peptide and dbSeqeunce data is loaded, this method creates
	 * the intermediate structure.
	 * 
	 * @return
	 */
	public IntermediateStructure buildIntermediateStructure() {
		if ((peptides.size() < 1) || (proteins.size() < 1)) {
			logger.error("no data to build the intermediate structure!");
			return null;
		}
		
		if (intermediateStructure != null) {
			logger.error("The intermediate structure was already created!");
			return null;
		}
		
        logger.info("creating intermediate structure with\n\t"
				+ getNrSpectrumIdentifications() + " spectrum identifications\n\t"
				+ getNrPeptides() + " peptides\n\t"
				+ getNrProteins() + " protein accessions");
        
		// first cluster the data
		List<Map<Comparable, Set<Comparable>>> clusterList = buildClusterList();
		
		// initialize the iterator
		clusterIterator = clusterList.listIterator();
		
		// initialize the intermediate structure
		intermediateStructure = new IntermediateStructure();
		
		// start the threads
		List<IntermediateStructureCreatorWorkerThread> threads;
		threads = new ArrayList<IntermediateStructureCreatorWorkerThread>(numberThreads);
		for (int i = 0; i < numberThreads; i++) {
			IntermediateStructureCreatorWorkerThread thread =
					new IntermediateStructureCreatorWorkerThread(i, this);
			threads.add(thread);
			
			thread.start();
		}
		
		// wait for the threads to finish
		for (IntermediateStructureCreatorWorkerThread thread : threads) {
			try {
				thread.join();
			} catch (InterruptedException e) {
				logger.error("thread got interrupted!", e);
				return null;
			}
		}
		
		logger.debug("intermediate structure contains "
				+ intermediateStructure.getNrClusters() + " clusters and "
				+ intermediateStructure.getNrGroups() + " groups");
		
		return intermediateStructure;
	}
	
	
	/**
	 * Creates mappings from peptide IDs to protein IDs, which are
	 * disjoint.
	 */
	private List<Map<Comparable, Set<Comparable>>> buildClusterList() {
		
		logger.info("start sorting clusters");
		
		// disjoint list of mappings from peptide IDs to protein IDs
		List<Map<Comparable, Set<Comparable>>> clusteredPepEntriesMap =
				new ArrayList<Map<Comparable,Set<Comparable>>>();
		
		Set<Comparable> peptidesDone = new HashSet<Comparable>(peptides.size());
		Set<Comparable> proteinsDone = new HashSet<Comparable>(proteins.size());
		
		for (Map.Entry<Comparable, Set<Comparable>> entryToPepsIt : proteinsToPeptidesMapping.entrySet()) {
			if (!proteinsDone.contains(entryToPepsIt.getKey())) {
				// this accession is not yet clustered, so start a new cluster
				// and insert all the "connected" peptides and accessions
				Map<Comparable, Set<Comparable>> pepEntriesMapCluster =
						createCluster(entryToPepsIt.getKey(), peptidesDone, proteinsDone);
				
				if (pepEntriesMapCluster != null) {
					clusteredPepEntriesMap.add(pepEntriesMapCluster);
				} else {
					logger.error("cluster could not be created!");
				}
				
				proteinsDone.add(entryToPepsIt.getKey());
			}
		}
		
		// the maps are no longer needed
		proteinsToPeptidesMapping = null;
		peptidesToProteinsMapping = null;
		
		logger.info(clusteredPepEntriesMap.size() + " sorted clusters");
		return clusteredPepEntriesMap;
	}
	
	
	/**
	 * Inserts the cluster of the given accession into the peptide accession
	 * map cluster.
	 * <p>
	 * This method should only be called by
	 * {@link IntermediateStructureCreator#buildClusterList()}
	 */
	private Map<Comparable, Set<Comparable>> createCluster(Comparable proteinID,
			Set<Comparable> peptidesDone, Set<Comparable> proteinsDone) {
		
		Set<Comparable> clusterProteins = new HashSet<Comparable>();
		Set<Comparable> clusterPeptides = new HashSet<Comparable>();
		
		boolean newProteins = true; // the given dbSequence is new 
		boolean newPeptides = false;
		
		// initialize the cluster's peptides with the peptides of the given dbSequence
		for (Comparable pepID : proteinsToPeptidesMapping.get(proteinID)) {
			newPeptides |= clusterPeptides.add(pepID);
		}
		
		// put the entries and peptides into the cluster and remove them from the original maps
		while (newProteins || newPeptides) {
			// repeat as long, as we get more accessions or peptides
			newProteins = false;
			newPeptides = false;
			
			// get proteins for peptides, which are new since the last loop
			for (Comparable newPeptideID : clusterPeptides) {
				if (!peptidesDone.contains(newPeptideID)) {
					for (Comparable newProteinID : peptidesToProteinsMapping.get(newPeptideID)) {
						newProteins |= clusterProteins.add(newProteinID);
					}
					peptidesDone.add(newPeptideID);
				}
			}
			
			// get peptides for proteins, which are new since the last loop
			for (Comparable newProteinID : clusterProteins) {
				if (!proteinsDone.contains(newProteinID)) {
					for (Comparable newPeptideID : proteinsToPeptidesMapping.get(newProteinID)) {
						newPeptides |= clusterPeptides.add(newPeptideID);
					}
					proteinsDone.add(newProteinID);
				}
			}
		}
		clusterProteins = null;
		
		// now we have the whole cluster, so put it into the pepAccMapCluster
		Map<Comparable, Set<Comparable>> peptidesToProteinsMapCluster =
				new HashMap<Comparable, Set<Comparable>>(clusterPeptides.size());
		
		for (Comparable peptideID : clusterPeptides) {
			peptidesToProteinsMapCluster.put(peptideID, peptidesToProteinsMapping.get(peptideID));
		}
		
		return peptidesToProteinsMapCluster;
	}
	
	
	/**
	 * Returns the next cluster in the clustered mapping of peptides to
	 * DBSequences.
	 * 
	 * @return
	 */
	protected synchronized Map<Comparable, Set<Comparable>> getNextCluster() {
		synchronized (clusterIterator) {
			if (clusterIterator != null) {
				if (clusterIterator.hasNext()) {
					return clusterIterator.next();
				} else {
					return null;
				}
			} else {
				logger.error("The cluster iterator is not yet initialized!");
				return null;
			}
		}
	}
	
	
	/**
	 * Adds the given cluster information into the intermediate structure.
	 * 
	 * @param cluster
	 */
	protected synchronized void addCluster(Collection<IntermediateGroup> cluster) {
		synchronized (intermediateStructure) {
			intermediateStructure.addCluster(cluster);
		}
	}
}
