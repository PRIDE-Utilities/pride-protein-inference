package uk.ac.ebi.pride.utilities.pia.modeller;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.jmztab.model.MZTabFile;
import uk.ac.ebi.pride.jmztab.utils.MZTabFileConverter;
import uk.ac.ebi.pride.utilities.data.controller.impl.ControllerImpl.MzIdentMLControllerImpl;
import uk.ac.ebi.pride.utilities.data.exporters.AbstractMzTabConverter;
import uk.ac.ebi.pride.utilities.data.exporters.HQMzIdentMLMzTabConverter;
import uk.ac.ebi.pride.utilities.data.exporters.MzIdentMLMzTabConverter;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateProtein;
import uk.ac.ebi.pride.utilities.pia.intermediate.prideimpl.PrideIntermediateProtein;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
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

import static org.junit.Assert.*;

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

            Map<Comparable, List<Comparable>> proteinPeptideMap = null;

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
}