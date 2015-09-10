package uk.ac.ebi.pride.utilities.pia.modeller.scores.protein;

import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptide;
import uk.ac.ebi.pride.utilities.pia.modeller.protein.inference.InferenceProteinGroup;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.ScoringItemType;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide.PeptideScoring;


/**
 * Calculates the score of the protein by adding up the scores of peptides.
 * 
 * @author julian
 *
 */
public class ProteinScoringAdditive extends ProteinScoring {
	
	public ProteinScoringAdditive(boolean countSpectrumOnce, PeptideScoring peptideScoring) {
		super(countSpectrumOnce, peptideScoring);
	}
	
	
	@Override
	public Double calculateProteinScore(InferenceProteinGroup proteinGroup) {
		Double proteinScore = Double.NaN;
		List<IntermediatePeptide> scorePeptides =
				new ArrayList<IntermediatePeptide>(proteinGroup.getPeptides());
		
		proteinGroup.removeAllScoringInformation();
		
//		if (countSpectrumOnce) {
//			// filter the peptides so that a spectrum counts only once
//			/*
//			Set<Comparable> scoringPSMsIDs = new HashSet<Comparable>();
//
//			Collections.sort(scorePeptides, new IntermediatePeptideComparator(peptideScoring));
//
//			ListIterator<IntermediatePeptide> pepIterator = scorePeptides.listIterator();
//			while (pepIterator.hasNext()) {
//				IntermediatePeptide peptide = pepIterator.next();
//
//				if ((peptide.getScore() != null) &&
//						!peptide.getScore().equals(Double.NaN)) {
//					boolean pepCounts = true;
//
//					List<IntermediatePeptideSpectrumMatch> pepsPSMs =
//							new ArrayList<IntermediatePeptideSpectrumMatch>(peptide.getPSMsWithScoringType(ScoringItemType.FULL_SCORING));
//					pepsPSMs.addAll(peptide.getPSMsWithScoringType(ScoringItemType.SHARED_SCORING));
//					for (IntermediatePeptideSpectrumMatch psm : pepsPSMs) {
//
//
//					}
//
//					if (pepCounts) {
//						// this peptide may be counted
//						for (IntermediatePeptideSpectrumMatch psm : pepsPSMs) {
//							scoringPSMsIDs.add(psm.getSpectrumIdentification().getSpectrum())
//
//						}
//
//
//
//					} else {
//
//					}
//
//
//				}
//
//			}
//			*/
//		}
		
		for (IntermediatePeptide peptide : scorePeptides) {
			if ((peptide.getScore() != null) &&
					!peptide.getScore().equals(Double.NaN)) {
				proteinGroup.setPeptidesScoringType(peptide, ScoringItemType.FULL_SCORING);
				// simply add up the scores of the peptides
				if (proteinScore.equals(Double.NaN)) {
					proteinScore = peptide.getScore();
				} else {
					proteinScore += peptide.getScore();
				}
			}
		}
		
		proteinGroup.setScore(proteinScore);
		return proteinScore;
	}
}