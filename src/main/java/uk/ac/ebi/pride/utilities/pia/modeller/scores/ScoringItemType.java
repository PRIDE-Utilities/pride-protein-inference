package uk.ac.ebi.pride.utilities.pia.modeller.scores;

public enum ScoringItemType {
	/** the item is fully used for the scoring */
	FULL_SCORING,
	
	/**
	 * the item is not used for scoring, because it shares another (sibling)
	 * item's score, which was selected for full scoring
	 */
	SHARED_SCORING,
	
	/** the item is not used for the scoring (e.g. due to filtering or another sibling is better suited) */
	NOT_SCORING,
	;
}
