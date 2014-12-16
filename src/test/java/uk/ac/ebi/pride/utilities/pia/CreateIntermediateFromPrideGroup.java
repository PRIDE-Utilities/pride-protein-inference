package uk.ac.ebi.pride.utilities.pia;

import static org.junit.Assert.*;

import java.io.File;
import java.net.URL;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.pride.utilities.data.controller.DataAccessController;
import uk.ac.ebi.pride.utilities.data.controller.impl.ControllerImpl.MzIdentMLControllerImpl;
import uk.ac.ebi.pride.utilities.pia.intermediate.DataImportController;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateGroup;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptide;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateStructure;
import uk.ac.ebi.pride.utilities.pia.intermediate.prideimpl.PrideImportController;
import uk.ac.ebi.pride.utilities.pia.modeller.PIAModeller;



public class CreateIntermediateFromPrideGroup {
	
	/** logger for this class */
	private static final Logger logger = LoggerFactory.getLogger(CreateIntermediateFromPrideGroup.class);
	
	/** the PRIDE access controller */
	private DataAccessController accessController = null;
	
	
	/**
	 * The expected results for each PAG.
	 * <p>
	 * The numbers for the groups are expected for each occurring intermediate
	 * group which has accessions. For proteins only the direct accessions,
	 * for peptides and PSMs all reachable.
	 * 
	 * @author julian
	 *
	 */
	enum ExpectedResults {
		PAG_hit_1 {
			@Override
			public int getNrProteins() {
				return 5;
			}
			
			@Override
			public int getNrPeptides() {
				return 80;
			}
			
			@Override
			public int getNrPSMs() {
				return 1090;
			}
			
			@Override
			public List<Integer> getNrProteinsInGroup() {
				return Arrays.asList(new Integer[]{1});
			}
			
			@Override
			public List<Integer> getNrPeptidesInGroup() {
				return Arrays.asList(new Integer[]{74, 7, 73, 72, 71});
			}
			
			@Override
			public List<Integer> getNrPSMsInGroup() {
				return Arrays.asList(new Integer[]{1077, 22, 1059, 1035, 1013});
			}
		},
		
		PAG_hit_2 {
			@Override
			public int getNrProteins() {
				return 2;
			}
			
			@Override
			public int getNrPeptides() {
				return 13;
			}
			
			@Override
			public int getNrPSMs() {
				return 28;
			}
			
			@Override
			public List<Integer> getNrProteinsInGroup() {
				return Arrays.asList(new Integer[]{1});
			}
			
			@Override
			public List<Integer> getNrPeptidesInGroup() {
				return Arrays.asList(new Integer[]{7, 10});
			}
			
			@Override
			public List<Integer> getNrPSMsInGroup() {
				return Arrays.asList(new Integer[]{17, 24});
			}
		},
		
		PAG_hit_3 {
			@Override
			public int getNrProteins() {
				return 30;
			}
			
			@Override
			public int getNrPeptides() {
				return 2;
			}
			
			@Override
			public int getNrPSMs() {
				return 2;
			}
			
			@Override
			public List<Integer> getNrProteinsInGroup() {
				return Arrays.asList(new Integer[]{4, 18, 8});
			}
			
			@Override
			public List<Integer> getNrPeptidesInGroup() {
				return Arrays.asList(new Integer[]{1, 2});
			}
			
			@Override
			public List<Integer> getNrPSMsInGroup() {
				return Arrays.asList(new Integer[]{1, 2});
			}
		},		
		
		PAG_hit_4 {
			@Override
			public int getNrProteins() {
				return 2;
			}
			
			@Override
			public int getNrPeptides() {
				return 1;
			}
			
			@Override
			public int getNrPSMs() {
				return 1;
			}
			
			@Override
			public List<Integer> getNrProteinsInGroup() {
				return Arrays.asList(new Integer[]{2});
			}
			
			@Override
			public List<Integer> getNrPeptidesInGroup() {
				return Arrays.asList(new Integer[]{1});
			}
			
			@Override
			public List<Integer> getNrPSMsInGroup() {
				return Arrays.asList(new Integer[]{1});
			}
		},		
		
		PAG_hit_5 {
			@Override
			public int getNrProteins() {
				return 1;
			}
			
			@Override
			public int getNrPeptides() {
				return 4;
			}
			
			@Override
			public int getNrPSMs() {
				return 8;
			}
			
			@Override
			public List<Integer> getNrProteinsInGroup() {
				return Arrays.asList(new Integer[]{1});
			}
			
			@Override
			public List<Integer> getNrPeptidesInGroup() {
				return Arrays.asList(new Integer[]{4});
			}
			
			@Override
			public List<Integer> getNrPSMsInGroup() {
				return Arrays.asList(new Integer[]{8});
			}
		},
		;
		
		
		public abstract int getNrProteins();
		
		public abstract int getNrPeptides();
		
		public abstract int getNrPSMs();
		
		public abstract List<Integer> getNrProteinsInGroup();
		
		public abstract List<Integer> getNrPeptidesInGroup();
		
		public abstract List<Integer> getNrPSMsInGroup();
	}
	
	@Before
	public void setUp() throws Exception {
		URL url = CreateIntermediateFromPrideGroup.class.getClassLoader().getResource("KPatel_111215_02.mzid");
		
		if (url == null) {
		    throw new IllegalStateException("no file for input found!");
		}
		
		File inputFile = new File(url.toURI());
		accessController = new MzIdentMLControllerImpl(inputFile);
		
		logger.debug("PAGs: " + accessController.getProteinAmbiguityGroupIds());
	}
	
	
	@After
	public void tearDown() throws Exception {
		if (accessController != null) {
			accessController.close();
		}
	}
	
	
	@Test
	public void createIntermediateGroupsFromAmbiguityGroups() {
		// ---------------------------------------------------------------------
		// creates intermediate structures for the protein ambiguity groups
		//
		
		for (ExpectedResults res : ExpectedResults.values()) {
			Comparable pag = res.toString();
			logger.debug("processing " + pag);
			
			PIAModeller piaModeller = new PIAModeller();
			
			piaModeller.addPrideControllerAsInput(accessController);
			DataImportController importController = piaModeller.getImportController(1);
			
			if (importController instanceof PrideImportController) {
				for (Comparable proteinId : accessController.getProteinAmbiguityGroupById(pag).getProteinIds()) {
					((PrideImportController) importController).addProteinsSpectrumIdentificationsToStructCreator(proteinId, null, null);
				}
			}
			
			IntermediateStructure intermediateStruct = piaModeller.buildIntermediateStructure();
			
			assertEquals("There should be only one cluster for " + pag, intermediateStruct.getNrClusters(), 1);
			
			assertEquals("Wrong number of proteins for " + pag, res.getNrProteins(), intermediateStruct.getNrProteins());
			assertEquals("Wrong number of peptides for " + pag, res.getNrPeptides(), intermediateStruct.getNrPeptides());
			assertEquals("Wrong number of PSMs for " + pag, res.getNrPSMs(), intermediateStruct.getAllIntermediatePSMs().size());
			
			for (IntermediateGroup group : intermediateStruct.getClusters().values().iterator().next()) {
				if ((group.getProteins() != null) && (group.getProteins().size() > 0)) {
					Integer nr = group.getProteins().size();
					assertTrue("Wrong number of proteins for group in " + pag + ", there should be no group with " + nr + " proteins",
							res.getNrProteinsInGroup().contains(nr));
					
					Set<IntermediatePeptide> peptides = group.getAllPeptides(); 
					nr = peptides.size();
					assertTrue("Wrong number of peptides for group in " + pag + ", there should be no group with " + nr + " peptides",
							res.getNrPeptidesInGroup().contains(nr));
					
					nr = 0;
					for (IntermediatePeptide pep : peptides) {
						nr += pep.getAllPeptideSpectrumMatches().size();
					}
					assertTrue("Wrong number of PSMs for group in " + pag + ", there should be no group with " + nr + " PSMs",
							res.getNrPSMsInGroup().contains(nr));
				}
			}
			
			piaModeller.close();
		}
	}
}
