package uk.ac.ebi.pride.utilities.pia.modeller.scores.protein;

import uk.ac.ebi.pride.utilities.pia.modeller.protein.inference.InferenceProteinGroup;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide.PeptideScoring;


public abstract class ProteinScoring {
	
	/** use a spectrum only once for protein scoring */
	protected boolean countSpectrumOnce;
	
	/** the peptide scoring used */
	protected PeptideScoring peptideScoring;
	
	
	/**
	 * 
	 * @param countSpectrumOnce if set to true and a spectrum is in the scoring
	 * set of more than one peptide, only the peptide with the highest score is
	 * used for full scoring. This resolves some conflicts with PSMs for
	 * different modifications from the same spectra
	 */
	public ProteinScoring(boolean countSpectrumOnce, PeptideScoring peptideScoring) {
		
		if (countSpectrumOnce) {
			System.err.println("Spectrum counting only once does not work yet!");
		}
		
		this.countSpectrumOnce = countSpectrumOnce;
		this.peptideScoring = peptideScoring;
	}
	
	
	/**
	 * Calculates the protein score for the inference group.
	 * 
	 * @param intermediatePeptide
	 */
	public abstract Double calculateProteinScore(InferenceProteinGroup proteinGroup);
}
