package uk.ac.ebi.pride.utilities.pia.modeller.psm;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.CvScore;


/**
 * This class handles the PSMs, originating from different importers (i.e.
 * input files).
 * 
 * @author julian
 *
 */
public class PSMModeller {
	
	/** logger for this class */
	private static final Logger logger = Logger.getLogger(PSMModeller.class);
	
	/** maps from the file (resp. controller) IDs to the corresponding PSMs */
	private Map<Integer, List<IntermediatePeptideSpectrumMatch>> filePSMs;
	
	/** maps from the fileIDs to all score accessions */
	private Map<Integer, Set<String>> fileScoreAccessions;
	
	/** mapping from the fileIDs to the score accessions used for FDR calculation */
	private Map<Integer, String> fileFdrScoreAccessions;
	
	
	
	public PSMModeller(Integer nrFiles) {
		filePSMs =
				new HashMap<Integer, List<IntermediatePeptideSpectrumMatch>>(nrFiles + 1);
		fileScoreAccessions = new HashMap<Integer, Set<String>>(nrFiles + 1);
		fileFdrScoreAccessions = new HashMap<Integer, String>(nrFiles + 1);
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
			fileScoreAccessions.put(fileID, new HashSet<String>(5));
		}
		
		/*
		TODO: 
		hier muessen alle scores vom PSM auch zu fileScoreAccessions geadded werden fileScoreAccessions.get(fileID);
		
		spectrumIdentification kann nicht für alle PSMs genommen werden... hat nur feste werte für scores!!!
		-> eigene klasse mit alle werten erstellen und SpectrumIdentification nur für entsprechende controller nehmen
		
		psm.getSpectrumIdentification().getScore().getCvTermReferenceWithValues()
		*/
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
		// TODO: wenn noch keine gesetzt, setze die jetzt automatisch auf (erste main) score
		
		return fileFdrScoreAccessions.get(fileID);
	}
	
	
	
	public String getFilesMainScoreAccession(Integer fileID) {
		if (fileScoreAccessions.get(fileID) == null) {
			return null;
		}
		
		String accession = null;
		
		for (String scoreAcc : fileScoreAccessions.get(fileID)) {
			CvScore cvScore = CvScore.getCvRefByAccession(scoreAcc);
			if ((cvScore != null) && cvScore.getMainScore()) {
				return cvScore.getAccession();
			}
			
			if (accession == null) {
				accession = scoreAcc;
			}
		}
		
		return accession;
	}
	
	
	public void calculateAllFDR() {
		
	}
	
	
	public void calculateFDR(Integer fileID) {
		
	}
}
