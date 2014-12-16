package uk.ac.ebi.pride.utilities.pia;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateProtein;
import uk.ac.ebi.pride.utilities.pia.modeller.PIAModeller;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterComparator;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.protein.ProteinAccessionFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.psm.PSMScoreFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.protein.inference.InferenceProteinGroup;
import uk.ac.ebi.pride.utilities.pia.modeller.protein.inference.OccamsRazorInference;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide.PeptideScoring;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide.PeptideScoringUseBestPSM;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.protein.ProteinScoring;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.protein.ProteinScoringAdditive;
import uk.ac.ebi.pride.utilities.term.CvTermReference;


public class ProteinInferenceTestPrideXML {
	
	/** logger for this class */
	private static final Logger logger = LoggerFactory.getLogger(ProteinInferenceTestPrideXML.class);
	
	
	private File inputFile;
	
	@Before
	public void setUp() throws Exception {
		inputFile = new File("/mnt/data/uniNOBACKUP/PrideInspector/manuscript-files/PXD000543/PRIDE_Exp_Complete_Ac_31941.xml");
	}
	
	@After
	public void tearDown() throws Exception {
		
	}
	
	@Test
	public void testWithoutImportfilter() throws Exception {
		PIAModeller piaModeller = new PIAModeller();
		
		// MS:1002355 		PSM-level FDRScore
		// MS:1002053 		MS-GF E-value
		// MS:1001171 		Mascot:score
		// PRIDE:0000069	Mascot:score
		//String scoreAcc = piaModeller.getPSMModeller().getFilesMainScoreAccession(1);
		String scoreAcc = "PRIDE:0000069";
		
		// add the input file to the files
		piaModeller.addFileAndImportSpectra(inputFile.getAbsolutePath(), null, scoreAcc);
		
		// first create the intermediate structure from the data given by the controller
		piaModeller.buildIntermediateStructure();
		
		// calculate FDR
		ProteinAccessionFilter decoyFilter = new ProteinAccessionFilter(FilterComparator.regex, "REV_.*", false);
		piaModeller.getPSMModeller().setDecoyFilter(decoyFilter);
		/*
		for (Integer id : piaModeller.getImportControllerIds()) {
			piaModeller.getPSMModeller().setFdrScoreAccession(id, 
					piaModeller.getPSMModeller().getFilesMainScoreAccession(id));
		}
		
		piaModeller.getPSMModeller().calculateAllFDR();
		*/
		
		
		// perform the protein inferences
		PeptideScoring pepScoring = new PeptideScoringUseBestPSM(scoreAcc, false);
		//ProteinScoring protScoring = new ProteinScoringMultiplicative(false, pepScoring);
		ProteinScoring protScoring = new ProteinScoringAdditive(false, pepScoring);
		
		
		List<AbstractFilter> filters = new ArrayList<AbstractFilter>();
		
		for (Double scoreLevel : new Double[]{0.0, 10.0/*, 20.0, 30.0, 40.0*/}) {
			
			logger.info("Processing score level " + scoreLevel);
			filters.clear();
			
			//filters.add(new PSMScoreFilter(FilterComparator.less_equal, scoreLevel, false, "MS:1002355", false));
			filters.add(new PSMScoreFilter(FilterComparator.greater_equal, scoreLevel, false, "PRIDE:0000069", false));
			
			
			piaModeller.getProteinModeller().infereProteins(
					pepScoring, protScoring, OccamsRazorInference.class, filters, false);
			/*
			piaModeller.getProteinModeller().infereProteins(
					pepScoring, protScoring, ReportAllInference.class, filters, false);
			*/
			// writing proteins to file
	        File outfile = new File(inputFile.getAbsolutePath() + ".mascot_score" + scoreLevel + ".txt");
	        BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(outfile));
	        
	        bufferedWriter.append("PAG\tPAG_score\t#peps\taccessions");
	        bufferedWriter.newLine();
			for (InferenceProteinGroup pag : piaModeller.getProteinModeller().getInferredProteins()) {
				// very simple export...
				StringBuilder protSB = new StringBuilder();
				for (IntermediateProtein prot : pag.getProteins()) {
					protSB.append(prot.getAccession());
					protSB.append(",");
				}
				
				bufferedWriter.append("\"" + pag.getID() + "\"\t\"" +
						pag.getScore() + "\"\t\"" +
						pag.getPeptides().size() + "\"\t\"" + 
						protSB + "\"");
				bufferedWriter.newLine();
			}
			bufferedWriter.close();
		}
		
		// close the modeller in the end
		piaModeller.close();
	}
}
