package uk.ac.ebi.pride.utilities.pia.modeller.scores.protein;

import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptide;
import uk.ac.ebi.pride.utilities.pia.modeller.protein.inference.InferenceProteinGroup;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.ScoringItemType;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide.PeptideScoring;


/**
 * Calculates the score of the protein by multiplying up the scores of peptides
 * and returning the -log10 of the score (log for numeric reasons, minus to
 * get higherScoreBetter)
 * 
 * @author julian
 *
 */
public class ProteinScoringMultiplicative extends ProteinScoring {
	
	public ProteinScoringMultiplicative(boolean countSpectrumOnce, PeptideScoring peptideScoring) {
		super(countSpectrumOnce, peptideScoring);
	}
	
	
	@Override
	public Double calculateProteinScore(InferenceProteinGroup proteinGroup) {
		Double proteinScore = Double.NaN;
		List<IntermediatePeptide> scorePeptides =
				new ArrayList<IntermediatePeptide>(proteinGroup.getPeptides());
		
		proteinGroup.removeAllScoringInformation();
		
		for (IntermediatePeptide peptide : scorePeptides) {
			if ((peptide.getScore() != null) &&
					!peptide.getScore().equals(Double.NaN)) {
				proteinGroup.setPeptidesScoringType(peptide, ScoringItemType.FULL_SCORING);
				
				// add up the -log10 of the scores
				if (proteinScore.equals(Double.NaN)) {
					proteinScore = -Math.log10(peptide.getScore());
				} else {
					proteinScore -= Math.log10(peptide.getScore());
				}
			}
		}
		
		proteinGroup.setScore(proteinScore);
		return proteinScore;
	}
}