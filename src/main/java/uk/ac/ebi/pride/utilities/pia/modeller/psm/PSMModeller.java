package uk.ac.ebi.pride.utilities.pia.modeller.psm;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.pia.modeller.fdr.FDRUtilities;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.protein.ProteinAccessionFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.CvScore;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.psm.IntermediatePSMComparator;



/**
 * This class handles the PSMs, originating from different importers (i.e.
 * input files).
 * 
 * @author julian
 *
 */
public class PSMModeller {
	
	/** logger for this class */
	private static final Logger logger =  LoggerFactory.getLogger(PSMModeller.class);
	
	
	/** maps from the file (resp. controller) IDs to the corresponding PSMs */
	private Map<Integer, List<IntermediatePeptideSpectrumMatch>> filePSMs;
	
	/** maps from the fileIDs to all score accessions */
	private Map<Integer, Set<String>> fileScoreAccessions;
	
	/** mapping from the fileIDs to the score accessions used for FDR calculation */
	private Map<Integer, String> fileFdrScoreAccessions;
	
	/** holds a filter to detect decoys. if no filter is given, the default decoy values are used for FDR calculation etc. */
	private ProteinAccessionFilter decoyFilter;
	
	/** whether to look for unknown CVs in the online OBO */
	private boolean oboLookup;
	
	
	public PSMModeller(Integer nrFiles, boolean oboLookup) {
		filePSMs = new HashMap<Integer, List<IntermediatePeptideSpectrumMatch>>(nrFiles + 1);
		fileScoreAccessions = new HashMap<Integer, Set<String>>(nrFiles + 1);
		fileFdrScoreAccessions = new HashMap<Integer, String>(nrFiles + 1);
		decoyFilter = null;
		
		this.oboLookup = oboLookup;
	}
	
	
	/**
	 * Adds the given PSM to the end of the list of PSMs for the file specified
	 * by fileID.
	 * 
	 * @param fileID
	 * @param psm
	 * @return
	 */
	public boolean addPSMforFile(Integer fileID, IntermediatePeptideSpectrumMatch psm) {
		if (!filePSMs.containsKey(fileID)) {
			filePSMs.put(fileID, new ArrayList<IntermediatePeptideSpectrumMatch>(10000));
		}
		
		return filePSMs.get(fileID).add(psm);
	}
	
	
	/**
	 * Returns the number of PSMs for the given file.
	 * 
	 * @param fileID
	 */
	public int getNrPSMs(Integer fileID) {
		if (!filePSMs.containsKey(fileID) || (filePSMs.get(fileID) == null)) {
			return 0;
		}
		
		return filePSMs.get(fileID).size();
	}
	
	
	/**
	 * Sets the given accession as the accession used for FDR calculation of the
	 * file given by fileID.
	 * 
	 * @param fileID
	 * @param accession
	 */
	public void setFdrScoreAccession(Integer fileID, String accession) {
		fileFdrScoreAccessions.put(fileID, accession);
	}
	
	
	/**
	 * Getter for the score set for FDR estimation of the file given by fileID
	 * @param fileID
	 * @return
	 */
	public String getFdrScoreAccession(Integer fileID) {
		if (fileFdrScoreAccessions.get(fileID) == null) {
			// if no score for FDR calculation is set yet, set the first main score
			fileFdrScoreAccessions.put(fileID, getFilesMainScoreAccession(fileID));
		}
		
		return fileFdrScoreAccessions.get(fileID);
	}
	
	
	/**
	 * Returns the main score of the file given by its ID.
	 * 
	 * @param fileID
	 * @return
	 */
	public String getFilesMainScoreAccession(Integer fileID) {
		if (!filePSMs.containsKey(fileID)) {
			return null;
		}
		
		String accession = null;
		for (String scoreAcc : getFilesScoreAccessions(fileID)) {
			CvScore cvScore = CvScore.getCvRefByAccession(scoreAcc);
			if ((cvScore != null) && cvScore.getIsMainScore()) {
				return cvScore.getAccession();
			}
			
			if (accession == null) {
				accession = scoreAcc;
			}
		}
		
		return accession;
	}
	
	
	/**
	 * Returns all possible score IDs of the file. This requires going through
	 * the PSMs at first call and might take some time.
	 * 
	 * @param fileID
	 * @return
	 */
	private Set<String> getFilesScoreAccessions(Integer fileID) {
		if (!filePSMs.containsKey(fileID)) {
			// this is an invalid file!
			return null;
		}
		
		if (fileScoreAccessions.get(fileID) == null) {
			// the score accessions are not yet indexed, do this no
			logger.info("getting the scores for file " + fileID);
			
			Set<String> accessions = new HashSet<String>();
			for (IntermediatePeptideSpectrumMatch psm : filePSMs.get(fileID)) {
				accessions.addAll(psm.getScoreAccessions());
			}
			fileScoreAccessions.put(fileID, accessions);
		}
		
		return fileScoreAccessions.get(fileID);
	}
	
	
	/**
	 * Returns the List of PSMs for the given file
	 * @param fileID
	 * @return
	 */
	public List<IntermediatePeptideSpectrumMatch> getFilesPSMs(Integer fileID) {
		return filePSMs.get(fileID);
	}
	
	
	/**
	 * Sets the filter for decoy identification.
	 */
	public void setDecoyFilter(ProteinAccessionFilter decoyFilter) {
		this.decoyFilter = decoyFilter;
	}
	
	
	/**
	 * Calculates the FDR of the PSMs for the file selected by its ID
	 */
	public void calculateFDR(Integer fileID) {
		List<IntermediatePeptideSpectrumMatch> psms = getFilesPSMs(fileID);
		String fdrScoreAccession = getFdrScoreAccession(fileID);
		
		FDRUtilities.markDecoys(psms, decoyFilter);
		
		Collections.sort(psms, new IntermediatePSMComparator(fdrScoreAccession, oboLookup));
		logger.info("PSMs sorted");
		
		FDRUtilities.calculateFDR(psms, fdrScoreAccession);
		logger.info("FDR calculated");
		
		FDRUtilities.calculateFDRScore(psms,  fdrScoreAccession, oboLookup);
		logger.info("FDR Score calculated");
	}
	
	
	/**
	 * Calculates the FDR of the PSMs for all files
	 */
	public void calculateAllFDR() {
		for (Integer fileID : filePSMs.keySet()) {
			calculateFDR(fileID);
		}
	}
}
