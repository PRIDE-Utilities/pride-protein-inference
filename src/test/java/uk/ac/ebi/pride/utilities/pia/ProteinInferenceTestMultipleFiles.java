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
import uk.ac.ebi.pride.utilities.pia.modeller.filter.peptide.PeptideNrPSMsFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.protein.ProteinNrPSMsFilter;
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
        int allowedThreads = 4;
        boolean considerModifications = false;
        String fdrScoreAccession = CvScore.PSI_MASCOT_SCORE.getAccession();
        String peptideScoreAccession = CvScore.PSI_MASCOT_SCORE.getAccession();
        boolean oboLookup = false;


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
        piaModeller.addFile("/mnt/data/uniNOBACKUP/PIA/testfiles/55merge/55merge_mascot_full.mzid");
        piaModeller.addFile("/mnt/data/uniNOBACKUP/PIA/testfiles/55merge/55merge_omssa.mzid");


        piaModeller.buildIntermediateStructure();


        // ---------------------------------------------------------------------
        // calculate FDR

        // sort the PSMs by score
		/*
		logger.info("sorting PSMs by score");
		List<IntermediatePeptideSpectrumMatch> psms =
				new ArrayList<IntermediatePeptideSpectrumMatch>(intermediateStructure.getAllIntermediatePSMs());
		logger.info("   obtained PSMs for sorting");

		Collections.sort(psms, new IntermediatePSMComparator(fdrScoreAccession, oboLookup));
		logger.info("   sorting done");

		// then calculate the FDR and FDRScore

		PSMAccessionsFilter decoyFilter =
				new PSMAccessionsFilter(FilterComparator.regex_only, "Rnd.*", false);
		int nrDecoys = FDRUtilities.markDecoys(psms, decoyFilter);
		logger.info("   decoys marked (" + nrDecoys + "/" + psms.size() + ")");

		FDRUtilities.calculateFDR(psms, fdrScoreAccession);
		logger.info("   fdr calculation done");

		FDRUtilities.calculateFDRScore(psms,  fdrScoreAccession, false);
		logger.info("   FDRScore calculation done");

		psms = null;	// this list is no longer needed
		*/

        // ---------------------------------------------------------------------
        // perform the protein inference
        //
		/*
		PeptideScoring pepScoring =
				new PeptideScoringUseBestPSM(peptideScoreAccession, oboLookup);
		ProteinScoring protScoring =
				new ProteinScoringAdditive(false, pepScoring);
		//		new ProteinScoringMultiplicative(false, pepScoring);

		AbstractProteinInference proteinInference =
				new OccamsRazorInference(intermediateStructure, pepScoring, protScoring, filters, allowedThreads);

		List<InferenceProteinGroup> inferenceGroups =
				proteinInference.calculateInference(considerModifications);

		logger.info("inferred groups: " + inferenceGroups.size());
        */

        // ---------------------------------------------------------------------
        // create the ProteinGroups from the inference groups
        //   -> these are the mzIdentML representatives of the groups (ProteinAmbiguityGroup)
        //
		/*
		List<ProteinGroup> proteinGroups =
				proteinInference.createProteinGroups(inferenceGroups);

		logger.info("Protein groups: " + proteinGroups.size());
		*/


        // ---------------------------------------------------------------------
        // print out some information
        //
		/*
		for (InferenceProteinGroup group : inferenceGroups) {

			StringBuilder groupText = new StringBuilder();

			for (IntermediateProtein protein : group.getProteins()) {
				if (groupText.length() > 0) {
					groupText.append(',');
				}
				groupText.append(protein.getAccession());
			}

			Set<String> subAccessions = new HashSet<String>();
			for (InferenceProteinGroup subGroup : group.getSubGroups()) {
				for (IntermediateProtein subProt : subGroup.getProteins()) {
					subAccessions.add(subProt.getAccession());
				}
			}
			for (String subAccession : subAccessions) {
				groupText.append(',');
				groupText.append('[');
				groupText.append(subAccession);
				groupText.append(']');
			}

			groupText.append(" (");
			groupText.append(group.getScore());
			groupText.append(')');
			///
			groupText.append('\n');

			for (IntermediatePeptide peptide : group.getPeptides()) {
				groupText.append('\t');
				groupText.append(peptide.getID() + ": " + peptide.getSequence() + " (" + peptide.getScore() + ", " + group.getPeptidesScoringType(peptide) + " )");


				for (IntermediatePeptideSpectrumMatch psm : peptide.getPeptideSpectrumMatches()) {
					groupText.append('\n');
					groupText.append('\t');
					groupText.append('\t');
					groupText.append(psm.getID() + ": " + psm.getScore(peptideScoreAccession) + " (" + peptide.getPSMsScoringType(psm) + ")");
				}

				groupText.append('\n');
			}
			///
			logger.info(groupText);
		}
		*/
		/*
        for (ProteinGroup group : proteinGroups) {

			StringBuilder accessions = new StringBuilder();
			Set<SpectrumIdentification> psmSet = new HashSet<SpectrumIdentification>();
			for (Protein protein : group.getProteinDetectionHypothesis()) {
				if (accessions.length() > 0) {
					accessions.append(',');
				}

				if (!protein.isPassThreshold()) {
					accessions.append('[');
					accessions.append(protein.getDbSequence().getAccession());
					accessions.append(']');
				} else {
					accessions.append(protein.getDbSequence().getAccession());
				}



				for (Peptide pep : protein.getPeptides()) {
					psmSet.add(pep.getSpectrumIdentification());
				}
			}

			StringBuilder psmSB = new StringBuilder();
			for (SpectrumIdentification psm : psmSet) {
				psmSB.append('\n');
				psmSB.append('\t');
				psmSB.append(psm.getId());
			}


			logger.info(accessions.toString() + psmSB.toString());
		}
		*/
    }
}