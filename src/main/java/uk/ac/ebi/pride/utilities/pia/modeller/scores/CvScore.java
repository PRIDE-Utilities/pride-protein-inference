package uk.ac.ebi.pride.utilities.pia.modeller.scores;

/**
 * Enumeration of known CV terms for scores
 * 
 * @author julian
 *
 */
public enum CvScore{
	
	PRIDE_OMSSA_E_VALUE("PRIDE", "PRIDE:0000185", "OMSSA E-value", false, false),
	PRIDE_OMSSA_P_VALUE("PRIDE", "PRIDE:0000186", "OMSSA P-value", false, true),
	PSI_OMSSA_E_VALUE("MS", "MS:1001328", "OMSSA:evalue", false, false),
	PSI_OMSSA_P_VALUE("MS", "MS:1001329", "OMSSA:pvalue", false, true),
	
	PRIDE_MASCOT_SCORE("PRIDE", "PRIDE:0000069", "Mascot Score", true, true),
	PRIDE_MASCOT_EXPECT_VALUE("PRIDE", "PRIDE:0000212", "Mascot expect value", false, false),
	PSI_MASCOT_SCORE("MS", "MS:1001171", "Mascot:score", true, true),
	PSI_MASCOT_EXPECT_VALUE("MS", "MS:1001172", "Mascot:expectation value", false, false),
	
	PRIDE_XTANDEM_HYPER_SCORE("PRIDE", "PRIDE:0000176", "X!Tandem Hyperscore", true, false),
	PRIDE_XTANDEM_EXPECTANCY_SCORE("PRIDE", "PRIDE:0000183", "X|Tandem expectancy score", false, true),
	PSI_XTANDEM_HYPERSCORE("MS", "MS:1001331", "X!Tandem:hyperscore", true, false),
	PSI_XTANDEM_EXPECTANCY_SCORE("MS", "MS:1001330", "X!Tandem:expect", false, true),
	
	PRIDE_SEQUEST_CN("PRIDE", "PRIDE:0000052", "Cn", true, false),
	PRIDE_SEQUEST_SCORE("PRIDE", "PRIDE:0000053", "SEQUEST SCORE", true, true),
	PRIDE_SEQUEST_DELTA_CN("PRIDE", "PRIDE:0000012", "Delta Cn", true, false),
	PSI_SEQUEST_CONSENSUS_SCORE("MS", "MS:1001163", "Sequest:consensus score", true, false),
	PSI_SEQUEST_DELTA_CN("MS", "MS:1001156", "Sequest:deltacn", true, false),
	PSI_SEQUEST_XCORR("MS", "MS:1001155", "Sequest:xcorr", true, true),
	
	PRIDE_PEPTIDE_PROPHET_DISCRIMINANT_SCORE("PRIDE", "PRIDE:0000138", "Discriminant score", true, false),
	PRIDE_PEPTIDE_PROPHET_PROBABILITY("PRIDE", "PRIDE:0000099", "PeptideProphet probability score", false, true),
	
	PSI_MYRIMATCH_MVH("MS", "MS:1001589", "MyriMatch:MVH", true, true),
	PSI_MYRIMATCH_NMATCHS("MS", "MS:1001121", "number of matched peaks", true, false),
	PSI_MYRIMATCH_NOMATCHS("MS", "MS:1001362", "number of unmatched peaks", false, false),
	PSI_MYRIMATCH_MZFIDELITY("MS", "MS:1001590", "MyriMatch:mzFidelity", true, false),
	
	//MS-GF
	PSI_MSGF_RAWSCORE("MS", "MS:1002049", "MS-GF raw score", true, false),
	PSI_MSGF_DENOVOSCORE("MS", "MS:1002050", "MS-GF de novo score", true, false),
	PSI_MSGF_SPECEVALUE("MS", "MS:1002052", "MS-GF spectral E-value", false, false),
	PSI_MSGF_EVALUE("MS", "MS:1002053", "MS-GF E-value", false, true),
	PSI_MSGF_QVALUE("MS", "MS:1002054", "MS-GF Q-value", false, false),
	PSI_MSGF_PEPQVALUE("MS", "MS:1002055", "MS-GF peptide-level Q-value", false, false),
	
	// Paragon
	//PSI_PARAGON_SCORE("MS", "MS:1001166", "Paragon:score", "MS:1001153"),
	
	// Phenyx
	PSI_PHENYX_SCORE("MS", "MS:1001390", "Phenyx:Score", true, true),
	
	// ProteinScape
	PSI_PROTEIN_EXTRACTOR_SCORE("MS", "MS:1001507", "ProteinExtractor:Score", true, false),
	PSI_PROTEINSCAPE_SEQUEST_METASCORE("MS", "MS:1001506", "ProteinScape:SequestMetaScore", true, false),
	
	// ProteinLynx
	//PSI_PROTEIN_LYNC_SCORE("MS", "MS:1001571", "ProteinLynx:Ladder Score", "MS:1001143"),
	
	// Sonar
	//PSI_SONAR_SCORE("MS", "MS:1001502", "Sonar:Score", "MS:1001143"),
	
	// percolator:score
	//PSI_PERCULATOR_SCORE("MS", "MS:1001492", "percolator:score", "MS:1001143"),
	
	
	PSI_PSM_LEVEL_LOCAL_FDR("MS", "MS:1002351", "PSM-level local FDR", false, false),
	PSI_PSM_LEVEL_Q_VALUE("MS", "MS:1002354", "PSM-level q-value", false, false),
	PSI_PSM_LEVEL_FDRSCORE("MS", "MS:1002355", "PSM-level FDRScore", false, false),
	PSI_PSM_LEVEL_COMBINED_FDRSCORE("MS", "MS:1002356", "PSM-level combined FDRScore", false, false),
	
	
	// PIA
	PSI_PIA_PROTEIN_SCORE("MS", "MS:1002394", "PIA:protein score", true, true),
	;
	
	
	private final String cvLabel;
	private final String accession;
	private final String name;
	private final boolean higherScoreBetter;
	private final boolean mainScore;
	
	private CvScore(String cvLabel, String accession, String name, boolean higherScoreBetter, boolean mainScore) {
		this.cvLabel = cvLabel;
		this.accession = accession;
		this.name = name;
		this.higherScoreBetter = higherScoreBetter;
		this.mainScore = mainScore;
	}
	
	
	public String getCvLabel() {
	    return cvLabel;
	}
	
	
	public String getAccession() {
	    return accession;
	}
	
	
	public String getName() {
		return name;
	}
	
	
	public boolean getHigherScoreBetter() {
		return higherScoreBetter;
	}
	
	
	public boolean getMainScore() {
		return mainScore;
	}
	
	
	/**
	 * Get Cv score by accession.
	 * @return the Cv score with the accession or null
	 */
	public static CvScore getCvRefByAccession(String accession) {
		for (CvScore cv : CvScore.values()) {
			if (cv.getAccession().equals(accession)) {
				return cv;
			}
		}
		
		return null;
	}
	
	
	/**
	 * Check whether the accession exists in the enum.
	 *
	 * @return boolean  true if exists
	 */
	public static boolean hasAccession(String accession) {
		for (CvScore cv : CvScore.values()) {
			if (cv.getAccession().equals(accession)) {
				return true;
			}
		}
		
		return false;
	}
}
