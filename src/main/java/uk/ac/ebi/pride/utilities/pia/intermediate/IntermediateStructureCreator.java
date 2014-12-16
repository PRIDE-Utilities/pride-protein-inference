package uk.ac.ebi.pride.utilities.pia.intermediate;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This class creates the intermediate structure needed for fast data access
 * during protein inference.
 * 
 * @author julian
 *
 */
public class IntermediateStructureCreator {
	
	/** the logger for this class */
	private static final Logger logger = LoggerFactory.getLogger(IntermediateStructureCreator.class);
	
	/** mapping from the peptide ID to intermediatePeptide
	 * TODO: we could decide in this class, whether a peptide is defined by the sequence only or also by the mods
	 **/
	private Map<Comparable, IntermediatePeptide> peptides;
	
	/** mapping from protein accession to the protein */
	private Map<String, IntermediateProtein> proteins;
	
	
	
	/** mapping from the protein accessions to connected peptides' IDs **/
	private Map<String, Set<Comparable>> proteinsToPeptidesMapping;
	
	/** mapping from the peptide IDs to connected proteins' accessions */
	private Map<Comparable, Set<String>> peptidesToProteinsMapping;
	
	
	
	/** iterates over the clustered list of peptide -> proteins mapping */
	private ListIterator<Map<Comparable, Set<String>>> clusterIterator;
	
	
	
	
	/** the created intermediate structure */
	private IntermediateStructure intermediateStructure;
	
	/** the maximal number of used threads */
	private int numberThreads;
	
	
	public IntermediateStructureCreator(int threads) {
		this.peptides = new HashMap<Comparable, IntermediatePeptide>();
		this.proteins = new HashMap<String, IntermediateProtein>();
		this.proteinsToPeptidesMapping =
				new HashMap<String, Set<Comparable>>();
		this.peptidesToProteinsMapping = new HashMap<Comparable, Set<String>>();
		
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
	 * given accessions
	 * @param proteinID
	 * @return
	 */
	public boolean proteinsContains(String accession) {
		return proteins.containsKey(accession);
	}
	
	
	/**
	 * adds the given protein to the proteins map
	 * 
	 * @param protein
	 * @return any previous protein with the same ID or null
	 */
	public IntermediateProtein addProtein(IntermediateProtein protein) {
		return proteins.put(protein.getAccession(), protein);
	}
	
	
	/**
	 * potentially adds further information to a protein, coming from e.g. a new
	 * file
	 * 
	 * @param proteinAccession the accession of the protein, which should get
	 * new information
	 * @param newProtein the protein, which may contain additional information
	 * @return any previous protein with the same ID or null
	 */
	public void addProteinInformation(String proteinAccession, IntermediateProtein newProtein) {
		// TODO: look if the newProtein contains additional information for the protein given by proteinAccession
		IntermediateProtein oldProtein = proteins.get(proteinAccession);
		
		if (oldProtein.getProteinSequence() == null) {
			if (newProtein.getProteinSequence() != null) {
				logger.warn("THIS IS NOT YET IMPLEMENTED!!! Don't use multiple files yet!");
			}
		} else if (newProtein.getProteinSequence() != null) {
			if (oldProtein.getProteinSequence() != newProtein.getProteinSequence()) {
				logger.warn("Protein with different sequences: " + oldProtein.getAccession() +
						", this is not supported, the first sequence is used" +
						"\n\t" + oldProtein.getProteinSequence() +
						"\n\t" + newProtein.getProteinSequence());
			}
		}
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
	public void addPeptideToProteinConnection(Comparable peptideID, String proteinAccession) {
		Set<String> protAccessions = peptidesToProteinsMapping.get(peptideID);
		if (protAccessions == null) {
			protAccessions = new HashSet<String>();
			peptidesToProteinsMapping.put(peptideID, protAccessions);
		}
		protAccessions.add(proteinAccession);
		
		Set<Comparable> pepIDs = proteinsToPeptidesMapping.get(proteinAccession);
		if (pepIDs == null) {
			pepIDs = new HashSet<Comparable>();
			proteinsToPeptidesMapping.put(proteinAccession, pepIDs);
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
		List<Map<Comparable, Set<String>>> clusterList = buildClusterList();
		
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
	 * Creates mappings from peptide IDs to protein accessions, which are
	 * disjoint.
	 */
	private List<Map<Comparable, Set<String>>> buildClusterList() {
		
		logger.info("start sorting clusters");
		
		// disjoint list of mappings from peptide IDs to protein IDs
		List<Map<Comparable, Set<String>>> clusteredPepEntriesMap =
				new ArrayList<Map<Comparable,Set<String>>>();
		
		Set<Comparable> peptidesDone = new HashSet<Comparable>(peptides.size());
		Set<String> proteinsDone = new HashSet<String>(proteins.size());
		
		for (Map.Entry<String, Set<Comparable>> entryToPepsIt : proteinsToPeptidesMapping.entrySet()) {
			if (!proteinsDone.contains(entryToPepsIt.getKey())) {
				// this accession is not yet clustered, so start a new cluster
				// and insert all the "connected" peptides and accessions
				Map<Comparable, Set<String>> pepEntriesMapCluster =
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
	private Map<Comparable, Set<String>> createCluster(String proteinAccession,
			Set<Comparable> peptidesDone, Set<String> proteinsDone) {
		
		Set<String> clusterProteins = new HashSet<String>();
		Set<Comparable> clusterPeptides = new HashSet<Comparable>();
		
		boolean newProteins = true; // the given dbSequence is new 
		boolean newPeptides = false;
		
		// initialize the cluster's peptides with the peptides of the given dbSequence
		for (Comparable pepID : proteinsToPeptidesMapping.get(proteinAccession)) {
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
					for (String newProteinID : peptidesToProteinsMapping.get(newPeptideID)) {
						newProteins |= clusterProteins.add(newProteinID);
					}
					peptidesDone.add(newPeptideID);
				}
			}
			
			// get peptides for proteins, which are new since the last loop
			for (String newProteinAccession : clusterProteins) {
				if (!proteinsDone.contains(newProteinAccession)) {
					for (Comparable newPeptideID : proteinsToPeptidesMapping.get(newProteinAccession)) {
						newPeptides |= clusterPeptides.add(newPeptideID);
					}
					proteinsDone.add(newProteinAccession);
				}
			}
		}
		clusterProteins = null;
		
		// now we have the whole cluster, so put it into the pepAccMapCluster
		Map<Comparable, Set<String>> peptidesToProteinsMapCluster =
				new HashMap<Comparable, Set<String>>(clusterPeptides.size());
		
		for (Comparable peptideID : clusterPeptides) {
			peptidesToProteinsMapCluster.put(peptideID, peptidesToProteinsMapping.get(peptideID));
		}
		
		return peptidesToProteinsMapCluster;
	}
	
	
	/**
	 * Returns the next cluster in the clustered mapping of peptides to
	 * proteins.
	 * 
	 * @return
	 */
	protected synchronized Map<Comparable, Set<String>> getNextCluster() {
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
