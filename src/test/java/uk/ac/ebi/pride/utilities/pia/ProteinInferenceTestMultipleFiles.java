package uk.ac.ebi.pride.utilities.pia;

import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import uk.ac.ebi.pride.utilities.pia.modeller.PIAModeller;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterComparator;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.psm.PSMDecoyFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.psm.PSMScoreFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.CvScore;


public class ProteinInferenceTestMultipleFiles {
	
	/** logger for this class */
	private static final Logger logger = Logger.getLogger(ProteinInferenceTestMultipleFiles.class);
	
	/** the PIA modeller */
	private PIAModeller piaModeller;
	
	
	@Before
	public void setUp() throws Exception {
		piaModeller = new PIAModeller();
	}
	
	@After
	public void tearDown() throws Exception {
		piaModeller.close();
	}
	
	@Test
	public void testProteinGroup() throws Exception {
		
		// ---------------------------------------------------------------------
		// some testing variables
		//

		boolean oboLookup = false;
		
		
		// ---------------------------------------------------------------------
		// set up some filters
		//
		List<AbstractFilter> filters = new ArrayList<AbstractFilter>();
		
		AbstractFilter filter =
				new PSMScoreFilter(FilterComparator.less_equal, 0.01, false, CvScore.PSI_PSM_LEVEL_FDRSCORE.getAccession(), oboLookup);
		filters.add(filter);
		
		filter = new PSMDecoyFilter(FilterComparator.equal, false, false);
		filters.add(filter);
		
		// ---------------------------------------------------------------------
		// first create the intermediate structure from the data given by the controller
		//
		piaModeller.addFile("/mnt/data/uniNOBACKUP/PIA/testfiles/55merge/55merge_mascot_full.mzid");
		piaModeller.addFile("/mnt/data/uniNOBACKUP/PIA/testfiles/55merge/55merge_omssa.mzid");
		
		
		piaModeller.buildIntermediateStructure();

	}
}
