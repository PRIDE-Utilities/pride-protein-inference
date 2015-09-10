package uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide;

import java.util.ArrayList;
import java.util.Iterator;
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
        Iterator<IntermediatePeptideSpectrumMatch> itPSM = intermediatePeptide.getPeptideSpectrumMatches().iterator();
		while(itPSM.hasNext()) {
            IntermediatePeptideSpectrumMatch psm = itPSM.next();
			Double score = psm.getScore(baseScoreAccession);
			
			if (bestScore.equals(Double.NaN) ||
					(psmComparator.compareValues(score, bestScore) <= 0)) {
				// the score may be null, if for some reason the PSM score could not be found
				if ((score != null) && (!bestScore.equals(score))) {
					bestScore = score;
					scoringPSM.clear();
				}
				scoringPSM.add(psm);
			}
		}
		
		if (scoringPSM.size() > 0) {
			boolean isFirst = true;
            Iterator<IntermediatePeptideSpectrumMatch> itMatch = scoringPSM.iterator();
			while(itMatch.hasNext()) {
                IntermediatePeptideSpectrumMatch psm = itMatch.next();
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
