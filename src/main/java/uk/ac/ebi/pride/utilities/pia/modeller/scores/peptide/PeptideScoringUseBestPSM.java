package uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide;

import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptide;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.ScoringItemType;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.psm.IntermediatePSMComparator;


/**
 * This class uses the PSM with the best score for scoring of a peptide.
 * 
 * @author julian
 *
 */
public class PeptideScoringUseBestPSM extends PeptideScoring {
	
	/** the used score comparator */
	protected IntermediatePSMComparator psmComparator;
	
	
	public PeptideScoringUseBestPSM(String scoreAccession, boolean oboLookup) {
		super(scoreAccession, oboLookup);
		
		psmComparator = new IntermediatePSMComparator(scoreAccession, oboLookup);
	}
	
	
	@Override
	public Double calculatePeptideScore(IntermediatePeptide intermediatePeptide) {
		Double bestScore = Double.NaN;
		List<IntermediatePeptideSpectrumMatch> scoringPSM = new ArrayList<IntermediatePeptideSpectrumMatch>();
		
		intermediatePeptide.removeAllScoringInformation();
		
		for (IntermediatePeptideSpectrumMatch psm : intermediatePeptide.getPeptideSpectrumMatches()) {
			Double score = psm.getScore(baseScoreAccession);
			
			if (bestScore.equals(Double.NaN) ||
					(psmComparator.compareValues(score, bestScore) <= 0)) {
				if (!bestScore.equals(score)) {
					bestScore = score;
					scoringPSM.clear();
				}
				scoringPSM.add(psm);
			}
		}
		
		if (scoringPSM.size() > 0) {
			boolean isFirst = true;
			for (IntermediatePeptideSpectrumMatch psm : scoringPSM) {
				// set just the first of the scoring PSMs to fully scoring
				intermediatePeptide.setPSMsScoringType(psm,
						isFirst ? ScoringItemType.FULL_SCORING : ScoringItemType.SHARED_SCORING);
				isFirst = false;
			}
			
			intermediatePeptide.setScore(bestScore);
			return bestScore;
		}
		
		intermediatePeptide.setScore(Double.NaN);
		return Double.NaN;
	}
	
	
}
