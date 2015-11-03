package uk.ac.ebi.pride.utilities.pia.modeller;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.jmztab.model.MZTabFile;
import uk.ac.ebi.pride.jmztab.utils.MZTabFileConverter;
import uk.ac.ebi.pride.utilities.data.controller.impl.ControllerImpl.MzIdentMLControllerImpl;
import uk.ac.ebi.pride.utilities.data.exporters.AbstractMzTabConverter;
import uk.ac.ebi.pride.utilities.data.exporters.HQMzIdentMLMzTabConverter;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateProtein;
import uk.ac.ebi.pride.utilities.pia.intermediate.prideimpl.PrideIntermediateProtein;
import uk.ac.ebi.pride.utilities.pia.modeller.protein.inference.InferenceProteinGroup;
import uk.ac.ebi.pride.utilities.pia.modeller.protein.inference.OccamsRazorInference;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.CvScore;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide.PeptideScoring;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide.PeptideScoringUseBestPSM;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.protein.ProteinScoring;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.protein.ProteinScoringAdditive;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.protein.ProteinScoringMultiplicative;
import uk.ac.ebi.pride.utilities.term.SearchEngineScoreCvTermReference;

import java.io.File;
import java.net.URL;
import java.util.*;

public class PIAModellerTest {

    MzIdentMLControllerImpl controller = null;

    @Before
    public void setUp() throws Exception {
        URL url = PIAModellerTest.class.getClassLoader().getResource("55merge_tandem.mzid");

        if (url == null) {
            throw new IllegalStateException("no file for input found!");
        }

        File inputFile = new File(url.toURI());
        controller = new MzIdentMLControllerImpl(inputFile);
    }

    @After
    public void tearDown() throws Exception {
        controller.close();
    }

    @Test
    public void runPIAModeller() throws  Exception{
        runDefaultProteinInference();
    }

    @Test
    public void runProteinInferenceToMzTab(){
        runDefaultProteinInference();
        AbstractMzTabConverter mzTabconverter = new HQMzIdentMLMzTabConverter(controller);
        MZTabFile mzTabFile = mzTabconverter.getMZTabFile();
        MZTabFileConverter checker = new MZTabFileConverter();
        checker.check(mzTabFile);
    }
    
    private void runDefaultProteinInference(){
        PIAModeller piaModeller = new PIAModeller();

        CvScore cvScore = null;

        String scoreAccession = null;
        // try to get the main-score
        for (SearchEngineScoreCvTermReference termRef : controller.getAvailablePeptideLevelScores()) {
            CvScore newCvScore;
            scoreAccession = termRef.getAccession();
            newCvScore = CvScore.getCvRefByAccession(termRef.getAccession());
            if ((newCvScore != null) && newCvScore.getIsMainScore()) {
                cvScore = newCvScore;
                scoreAccession = cvScore.getAccession();
                break;
            }
        }

        // add the input file to modeller and import data
        Integer controllerID = piaModeller.addPrideControllerAsInput(controller);

        piaModeller.importAllDataFromFile(controllerID);

        // first create the intermediate structure from the data given by the controller
        piaModeller.buildIntermediateStructure();

        PeptideScoring pepScoring = new PeptideScoringUseBestPSM(scoreAccession, false);
        ProteinScoring protScoring;
        if ((cvScore != null) && !cvScore.getHigherScoreBetter()) {
            protScoring = new ProteinScoringMultiplicative(false, pepScoring);
        } else {
            protScoring = new ProteinScoringAdditive(false, pepScoring);
        }

        // perform the protein inferences
        piaModeller.getProteinModeller().infereProteins(pepScoring, protScoring, OccamsRazorInference.class, null, false);

        // create the protein groups
        int nrGroups = piaModeller.getProteinModeller().getInferredProteins().size();
        Map<Comparable, Map<Comparable, List<Comparable>>> prideProteinGroupMapping = new HashMap<Comparable, Map<Comparable,List<Comparable>>>(nrGroups);

        for (InferenceProteinGroup piaGroup : piaModeller.getProteinModeller().getInferredProteins()) {

            Map<Comparable, List<Comparable>> proteinPeptideMap;

            Set<IntermediateProtein> proteinSet = new HashSet<IntermediateProtein>(piaGroup.getProteins());

            // include the subGroups
            for (InferenceProteinGroup subGroup : piaGroup.getSubGroups()) {
                proteinSet.addAll(subGroup.getProteins());
            }

            proteinPeptideMap = new HashMap<Comparable, List<Comparable>>(proteinSet.size());

            for (IntermediateProtein protein : proteinSet) {
                Comparable proteinID = ((PrideIntermediateProtein)protein).getPrideProteinID();
                proteinPeptideMap.put(proteinID, null);
            }
            prideProteinGroupMapping.put(piaGroup.getID(), proteinPeptideMap);
        }
        controller.setInferredProteinGroups(prideProteinGroupMapping);
    }
    
    
	@Test
	public void inferenceFromMultipleFiles() {
		// this test performs a protein inference on multiple MS-GF+ files
		
		// setting the main score for these data
		String scoreAccession = "MS:1002053";
		CvScore cvScore = CvScore.PSI_MSGF_EVALUE;
		
		MzIdentMLControllerImpl controller1 = null;
		MzIdentMLControllerImpl controller2 = null;
		MzIdentMLControllerImpl controller3 = null;
		
		try {
			URL url = PIAModellerTest.class.getClassLoader().getResource("PXD001428/F1_20150223_Agilent5_PG_incl_list_AspN_A[Node_05].mzid");
	        File inputFile = new File(url.toURI());
	        controller1 = new MzIdentMLControllerImpl(inputFile);
			
			url = PIAModellerTest.class.getClassLoader().getResource("PXD001428/F1_20150223_Agilent5_PG_incl_list_AspN_A[Node_08].mzid");
	        inputFile = new File(url.toURI());
	        controller2 = new MzIdentMLControllerImpl(inputFile);
	        
			url = PIAModellerTest.class.getClassLoader().getResource("PXD001428/F1_20150223_Agilent5_PG_incl_list_AspN_B[Node_05].mzid");
	        inputFile = new File(url.toURI());
	        controller3 = new MzIdentMLControllerImpl(inputFile);
		} catch (Exception e) {
			throw new IllegalStateException("Problem reading file!", e);
		}
		
		
		PIAModeller piaModeller = new PIAModeller();
		
		// add the input file to modeller and import data
		Integer controllerID1 = piaModeller.addPrideControllerAsInput(controller1);
		piaModeller.importAllDataFromFile(controllerID1);
		
		Integer controllerID2 = piaModeller.addPrideControllerAsInput(controller2);
		piaModeller.importAllDataFromFile(controllerID2);
		
		Integer controllerID3 = piaModeller.addPrideControllerAsInput(controller3);
		piaModeller.importAllDataFromFile(controllerID3);
		
		// first create the intermediate structure from the data given by the controller
		piaModeller.buildIntermediateStructure();
		
		// score the peptides
		PeptideScoring pepScoring = new PeptideScoringUseBestPSM(scoreAccession, false);
		ProteinScoring protScoring;
		if ((cvScore != null) && !cvScore.getHigherScoreBetter()) {
		    protScoring = new ProteinScoringMultiplicative(false, pepScoring);
		} else {
		    protScoring = new ProteinScoringAdditive(false, pepScoring);
		}
		
		// perform the protein inferences
		piaModeller.getProteinModeller().infereProteins(pepScoring, protScoring, OccamsRazorInference.class, null, false);
		
		// create the protein groups
		int nrGroups = piaModeller.getProteinModeller().getInferredProteins().size();
		Map<Comparable, Map<Comparable, List<Comparable>>> prideProteinGroupMapping = new HashMap<Comparable, Map<Comparable,List<Comparable>>>(nrGroups);
		
		for (InferenceProteinGroup piaGroup : piaModeller.getProteinModeller().getInferredProteins()) {
			Map<Comparable, List<Comparable>> proteinPeptideMap;
			
			Set<IntermediateProtein> proteinSet = new HashSet<IntermediateProtein>(piaGroup.getProteins());
			
			// include the subGroups
			for (InferenceProteinGroup subGroup : piaGroup.getSubGroups()) {
			    proteinSet.addAll(subGroup.getProteins());
			}
			
			proteinPeptideMap = new HashMap<Comparable, List<Comparable>>(proteinSet.size());
			
			for (IntermediateProtein protein : proteinSet) {
			    Comparable proteinID = ((PrideIntermediateProtein)protein).getPrideProteinID();
			    proteinPeptideMap.put(proteinID, null);
			}
			prideProteinGroupMapping.put(piaGroup.getID(), proteinPeptideMap);
		}
		
		// TODO: this needs another controller, everything below will not work
		/*
		controller1.setInferredProteinGroups(prideProteinGroupMapping);
		
		AbstractMzTabConverter mzTabconverter = new HQMzIdentMLMzTabConverter(controller1);
		MZTabFile mzTabFile = mzTabconverter.getMZTabFile();
		MZTabFileConverter checker = new MZTabFileConverter();
		checker.check(mzTabFile);
		*/
	}
}