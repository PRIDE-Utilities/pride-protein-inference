package uk.ac.ebi.pride.utilities.pia.modeller.scores.psm;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.ScoreUtilities;

/**
 * Compares the scores of {@link IntermediatePeptideSpectrumMatch}es
 * 
 * @author julian
 *
 */
public class IntermediatePSMComparator implements Comparator<IntermediatePeptideSpectrumMatch>{
	
	/** the accession of the compared score */
	private String scoreAccession;
	
	/** whether a higher score is better for the compared PSM score */
	private boolean higherScoreBetter;
	
	/** a score lookup map */
	private Map<Comparable, Double> scoreLookup;
	
	
	public IntermediatePSMComparator(String scoreAccession, boolean oboLookup) {
		this.scoreAccession = scoreAccession;
		this.higherScoreBetter = ScoreUtilities.isHigherScoreBetter(scoreAccession, oboLookup);
		this.scoreLookup = new HashMap<Comparable, Double>();
	}
	
	
	@Override
	public int compare(IntermediatePeptideSpectrumMatch o1,
			IntermediatePeptideSpectrumMatch o2) {
		Double score1;
		Double score2;
		
		if (o1 != null) {
			if (scoreLookup.containsKey(o1.getID())) {
				score1 = scoreLookup.get(o1.getID());
			} else {
				score1 = o1.getScore(scoreAccession);
				scoreLookup.put(o1.getID(), score1);
			}
		} else {
			score1 = null;
		}
		
		if (o2 != null) {
			if (scoreLookup.containsKey(o2.getID())) {
				score2 = scoreLookup.get(o2.getID());
			} else {
				score2 = o2.getScore(scoreAccession);
				scoreLookup.put(o2.getID(), score2);
			}
		} else {
			score2 = null;
		}
		
		return compareValues(score1, score2);
	}
	
	
	/**
	 * Compares the values of two PSM scores.
	 * 
	 * @param o1
	 * @param o2
	 * @return
	 */
	public int compareValues(Double score1, Double score2) {
		return ScoreUtilities.compareValues(score1, score2, higherScoreBetter);
	}
}
