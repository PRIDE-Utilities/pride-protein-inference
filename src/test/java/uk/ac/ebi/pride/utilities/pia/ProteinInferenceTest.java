package uk.ac.ebi.pride.utilities.pia;

import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;


import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.pride.utilities.pia.intermediate.DataImportController;

import uk.ac.ebi.pride.utilities.pia.modeller.PIAModeller;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterComparator;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.peptide.PeptideNrPSMsFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.protein.ProteinAccessionFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.protein.ProteinNrPSMsFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.psm.PSMDecoyFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.psm.PSMScoreFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.CvScore;

public class ProteinInferenceTest {
	
	/** logger for this class */
	private static final Logger logger = LoggerFactory.getLogger(ProteinInferenceTest.class);
	
	private File inputFile = null;
	
	private DataImportController importController = null;
	
	@Before
	public void setUp() throws Exception {
		URL url = ProteinInferenceTest.class.getClassLoader().getResource("report-proteins-report_all-55merge_mascot_full.mzid");

		if (url == null) {
		    throw new IllegalStateException("no file for input found!");
		}
		
		inputFile = new File(url.toURI());
	}
	
	@After
	public void tearDown() throws Exception {
		if (importController != null) {
			importController.close();
		}
	}
	
	@Test
	public void testProteinGroup() throws Exception {
		
		// ---------------------------------------------------------------------
		// some testing variables
		//
		int allowedThreads = 4;
		boolean considerModifications = false;
		boolean filterPSMsOnImport = false;
		String fdrScoreAccession = CvScore.PSI_MASCOT_SCORE.getAccession();
		String peptideScoreAccession = CvScore.PSI_MASCOT_SCORE.getAccession();
		boolean oboLookup = false;
		
		
		logger.info("using " + inputFile.getAbsolutePath());
		
		// ---------------------------------------------------------------------
		// set up some filters
		//
		List<AbstractFilter> filters = new ArrayList<AbstractFilter>();
		
		AbstractFilter filter =
		//		new PSMScoreFilter(FilterComparator.greater_equal, 20.0, false, CvScore.PSI_MASCOT_SCORE.getAccession(), oboLookup);
				new PSMScoreFilter(FilterComparator.less_equal, 0.01, false, CvScore.PSI_PSM_LEVEL_FDRSCORE.getAccession(), oboLookup);
		filters.add(filter);
		
		filter = new PSMDecoyFilter(FilterComparator.equal, false, false);
		filters.add(filter);
		
		filter = new PeptideNrPSMsFilter(FilterComparator.greater_equal, 2, false);
		//filters.add(filter);
		
		filter = new ProteinNrPSMsFilter(FilterComparator.greater_equal, 4, false);
		//filters.add(filter);
		
		
		// ---------------------------------------------------------------------
		// first create the intermediate structure from the data given by the controller
		//
		PIAModeller piaModeller = new PIAModeller();
		
        if (filterPSMsOnImport) {
        	piaModeller.addFileAndImportSpectra(inputFile.getAbsolutePath(), filters);
        } else {
        	piaModeller.addFileAndImportSpectra(inputFile.getAbsolutePath(), null);
        }
		
        piaModeller.buildIntermediateStructure();
		
        
		// ---------------------------------------------------------------------
		// calculate FDR
		ProteinAccessionFilter decoyFilter =
				new ProteinAccessionFilter(FilterComparator.regex, "Rnd.*", false);
		
		piaModeller.getPSMModeller().setDecoyFilter(decoyFilter);
		
		piaModeller.getPSMModeller().calculateAllFDR();
		
		piaModeller.close();
	}
}
