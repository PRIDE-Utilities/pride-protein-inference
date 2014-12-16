package uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide;

import java.util.Comparator;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptide;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.ScoreUtilities;



/**
 * Compares the scores of {@link IntermediatePeptide}s. These have to be
 * calculated beforehand.
 * 
 * @author julian
 *
 */
public class IntermediatePeptideComparator implements Comparator<IntermediatePeptide> {
	
	/** whether a higher score is better for the compared peptide score */
	private boolean higherScoreBetter;
	
	
	/**
	 * Creates the comparator with the given value for higherScoreBetter.
	 * 
	 * @param higherScoreBetter
	 */
	public IntermediatePeptideComparator(boolean higherScoreBetter) {
		this.higherScoreBetter = higherScoreBetter;
	}
	
	
	/**
	 * Creates the comparator and gets the value for higherScoreBetter from the
	 * given {@link PeptideScoring}.
	 * 
	 * @param higherScoreBetter
	 */
	public IntermediatePeptideComparator(PeptideScoring peptideScoring) {
		this.higherScoreBetter = peptideScoring.higherScoreBetter;
	}
	
	
	@Override
	public int compare(IntermediatePeptide o1, IntermediatePeptide o2) {
		Double score1;
		Double score2;
		
		if (o1 != null) {
			score1 = o1.getScore();
		} else {
			score1 = null;
		}
		
		if (o2 != null) {
			score2 = o2.getScore();
		} else {
			score2 = null;
		}
		
		return compareValues(score1, score2);
	}
	
	
	/**
	 * Compares the values of two peptide scores.
	 * 
	 * @param o1
	 * @param o2
	 * @return
	 */
	public int compareValues(Double score1, Double score2) {
		return ScoreUtilities.compareValues(score1, score2, higherScoreBetter);
	}
}
