package uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptide;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.CvScore;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.ScoreUtilities;


public abstract class PeptideScoring {
	
	/** accession of the base score */
	protected String baseScoreAccession;
	
	/** whether higher value of the base score is better than a lower value */
	protected boolean higherScoreBetter;
	
	
	/**
	 * Creates a new scoring object and determines (if appropriate), if a higher
	 * value is better for the given base-score. 
	 * 
	 * @param scoreAccession
	 * @param oboLookup
	 */
	public PeptideScoring(String scoreAccession, boolean oboLookup) {
		if (!CvScore.hasAccession(scoreAccession)) {
			if (!oboLookup ||
					(ScoreUtilities.findAccessionInObo(scoreAccession) == null)) {
				scoreAccession = null;
			}
		}
		
		this.baseScoreAccession = scoreAccession;
		higherScoreBetter = ScoreUtilities.isHigherScoreBetter(scoreAccession, oboLookup);
	}
	
	
	/**
	 * Calculates the peptide score for the intermediate peptide. For the
	 * calculation, the PSMs which are retrieved from the intermediate peptide,
	 * may be only the ones passing a set filter, if the peptide was filtered
	 * before.
	 * <p>
	 * The scoring PSMs are set in the intermediatePeptide, as well as the final
	 * score.
	 * 
	 * @param intermediatePeptide
	 */
	public abstract Double calculatePeptideScore(IntermediatePeptide intermediatePeptide);
}
