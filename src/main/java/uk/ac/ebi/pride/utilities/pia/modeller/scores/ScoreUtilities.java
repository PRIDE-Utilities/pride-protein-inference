package uk.ac.ebi.pride.utilities.pia.modeller.scores;

public class ScoreUtilities {
	/**
	 * Checks the PSI-MS obo and returns the name of the accession.
	 * 
	 * @param scoreAccession
	 * @return name of the accession or null, if not found
	 */
	public static String findAccessionInObo(String scoreAccession) {
		// TODO: implement!
		return null;
	}
	
	
	/**
	 * Returns for a given score accession, whether a higher score is better.
	 * <p>
	 * If the score is not known and oboLookup is true, the score is searched
	 * in the PSI-MS obo.
	 * <p>
	 * If nothing about the score is found, true will be returned.
	 * 
	 * @param scoreAccession
	 */
	public static boolean isHigherScoreBetter(String scoreAccession, boolean oboLookup) {
		
		CvScore cvScore = CvScore.getCvRefByAccession(scoreAccession);
		if (cvScore != null) {
			return cvScore.getHigherScoreBetter();
		}
		
//		if (oboLookup) {
//			// TODO: implement
//		}
		
		return false;
	}
	
	
	/**
	 * Compares the values of two values considered to be scores.
	 * <p>
	 * null and NaN is always considered to be worse than any other value.
	 * 
	 * @param score1
	 * @param score2
	 * @param higherScoreBetter whether a higher score value is considered to be better
	 * @return an int smaller 0, if score1 is considered better than score 2, 0
	 * if both are considered to be equal, bigger 0 if score2 is better than
	 * score1
	 */
	public static int compareValues(Double score1, Double score2, boolean higherScoreBetter) {
		if ((score1 == null) && 
				(score2 == null)) {
			// both PSMs don't have this score
			return 0;
		} else if ((score1 == null) && 
				(score2 != null)) {
			// psm1 does not have the score
			return 1;
		} else if ((score1 != null) && 
				(score2 == null)) {
			// psm2 does not have the score
			return -1;
		} else {
			// both scores are != null
			if (score1 != null ? score1.equals(Double.NaN) : false) {
				// NaN is always considered worse than anything (null likewise)
				return 1;
			} else if (score2 != null ? score2.equals(Double.NaN) : false) {
				return -1;
			} else {
				if (higherScoreBetter) {
					return -(score1 != null ? score1.compareTo(score2) : 0);
				} else {
					return score1 != null ? score1.compareTo(score2) : 0;
				}
			}
		}
	}
}
