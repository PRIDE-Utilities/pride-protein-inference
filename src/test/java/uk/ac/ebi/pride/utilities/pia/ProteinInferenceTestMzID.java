package uk.ac.ebi.pride.utilities.pia;

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
import uk.ac.ebi.pride.utilities.pia.modeller.protein.inference.ReportAllInference;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide.PeptideScoring;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide.PeptideScoringUseBestPSM;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.protein.ProteinScoring;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.protein.ProteinScoringMultiplicative;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

public class ProteinInferenceTestMzID {

    /** logger for this class */
    private static final Logger logger = LoggerFactory.getLogger(ProteinInferenceTestMzID.class);


    private File inputFile;

    @Before
    public void setUp() throws Exception {
        //inputFile = new File("/mnt/data/uniNOBACKUP/PrideInspector/manuscript-files/PXD000320/QC_12_03_1000ng_run09_15Feb13_Rogue_12-10-04_msgfplus.mzid");
        //inputFile = new File("/mnt/data/uniNOBACKUP/PrideInspector/manuscript-files/PXD000320/QC_Shew_LS_30ug_500min_2_msgfplus.mzid");

        inputFile = new File("/mnt/data/uniNOBACKUP/PrideInspector/manuscript-files/PXD000320/QC_Shew_12_02_Run-01_18Sep12_Eagle_12-06-09_msgfplus.mzid");
    }

    @After
    public void tearDown() throws Exception {

    }

    @Test
    public void testWithoutImportfilter() throws Exception {
        PIAModeller piaModeller = new PIAModeller();

        // MS:1002355 	PSM-level FDRScore
        // MS:1002053 	MS-GF E-value
        // MS:1001171 	Mascot:score
        //String scoreAcc = piaModeller.getPSMModeller().getFilesMainScoreAccession(1);
        String scoreAcc = "MS:1002053";


        // add the input file to the files
        piaModeller.addFileAndImportSpectra(inputFile.getAbsolutePath(), null, scoreAcc);

        // first create the intermediate structure from the data given by the controller
        piaModeller.buildIntermediateStructure();

        // calculate FDR
        ProteinAccessionFilter decoyFilter = new ProteinAccessionFilter(FilterComparator.regex, "XXX_.*", false);
        piaModeller.getPSMModeller().setDecoyFilter(decoyFilter);

        // perform the protein inferences
        PeptideScoring pepScoring = new PeptideScoringUseBestPSM(scoreAcc, false);
        ProteinScoring protScoring = new ProteinScoringMultiplicative(false, pepScoring);


        List<AbstractFilter> filters = new ArrayList<AbstractFilter>();
        // report-all
        for (Double scoreLevel : new Double[]{0.01, 0.03, 0.05}) {

            logger.info("Processing score level " + scoreLevel);
            filters.clear();

            filters.add(new PSMScoreFilter(FilterComparator.less_equal, scoreLevel, false, "MS:1002053", false));

            piaModeller.getProteinModeller().infereProteins(
                    pepScoring, protScoring, ReportAllInference.class, filters, false);

            // writing proteins to file
            File outfile = new File(inputFile.getAbsolutePath() + ".report_all.msgf_evalue" + scoreLevel + ".txt");
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

        // occams-razor
        for (Double scoreLevel : new Double[]{0.01, 0.03, 0.05}) {

            logger.info("Processing score level " + scoreLevel);
            filters.clear();

            filters.add(new PSMScoreFilter(FilterComparator.less_equal, scoreLevel, false, "MS:1002053", false));

            piaModeller.getProteinModeller().infereProteins(
                    pepScoring, protScoring, OccamsRazorInference.class, filters, false);

            // writing proteins to file
            File outfile = new File(inputFile.getAbsolutePath() + ".occams.msgf_evalue" + scoreLevel + ".txt");
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

        // occams-razor and FDR-score
        piaModeller.getPSMModeller().setFdrScoreAccession(1, "MS:1002053");
        piaModeller.getPSMModeller().calculateAllFDR();

        scoreAcc = "MS:1002355";
        pepScoring = new PeptideScoringUseBestPSM(scoreAcc, false);
        protScoring = new ProteinScoringMultiplicative(false, pepScoring);

        for (Double scoreLevel : new Double[]{0.01, 0.03, 0.05}) {

            logger.info("Processing score level " + scoreLevel);
            filters.clear();

            filters.add(new PSMScoreFilter(FilterComparator.less_equal, scoreLevel, false, "MS:1002355", false));

            piaModeller.getProteinModeller().infereProteins(
                    pepScoring, protScoring, OccamsRazorInference.class, filters, false);

            // writing proteins to file
            File outfile = new File(inputFile.getAbsolutePath() + ".occams.fdrscore" + scoreLevel + ".txt");
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