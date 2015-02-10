package uk.ac.ebi.pride.utilities.pia.intermediate.prideimpl;

import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.jmzidml.model.mzidml.AbstractParam;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.pride.utilities.data.controller.DataAccessController;
import uk.ac.ebi.pride.utilities.data.core.SpectrumIdentification;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.CvScore;
import uk.ac.ebi.pride.utilities.term.CvTermReference;


/**
 * Representation of a peptide spectrum match in the intermediate structure.
 * 
 * @author julian
 *
 */
public class PrideIntermediatePeptideSpectrumMatch extends IntermediatePeptideSpectrumMatch {
	
	/** a unique ID (cached on the first accession) */
	private String id;
	
	/** the used PRIDE dataAccessController */
	private DataAccessController controller;
	
	/** the protein ID for accession by a PRIDE dataAccessController */
	private Comparable proteinID;
	
	/** the peptide ID for accession by a PRIDE dataAccessController */
	private Comparable peptideID;
	
	/** the accession of a cached score */
	private String cachedScoreAccession;
	
	/** the value of a cached score */
	private Double cachedScore;
	
	
	public PrideIntermediatePeptideSpectrumMatch(DataAccessController controller,
			Comparable proteinID, Comparable peptideID) {
		this(controller, proteinID, peptideID, null);
	}
	
	
	/**
	 * creates a new PRIDE intermediate PSM and caches the score given by the
	 * accession (this requires some mor memory, but is much faster)
	 * 
	 * @param controller
	 * @param proteinID
	 * @param peptideID
	 * @param cacheScoreAccession
	 */
	public PrideIntermediatePeptideSpectrumMatch(DataAccessController controller,
			Comparable proteinID, Comparable peptideID, String cacheScoreAccession) {
		super();
		
		this.id = null;
		this.controller = controller;
		this.proteinID = proteinID;
		this.peptideID = peptideID;
		
		// get the score and cache it
		if (cacheScoreAccession != null) {
			Double score = getScore(cacheScoreAccession);
			if (score != null) {
				this.cachedScore = score;
				this.cachedScoreAccession = cacheScoreAccession;
			}
		} else {
			this.cachedScore = null;
			this.cachedScoreAccession = null;
		}
	}
	
	
	@Override
	public String getID() {
		if (id == null) {
			id = getControllerID() + ":" + getSpectrumIdentification().getId();
		}
		return id;
	}
	
	
	@Override
	public String getControllerID() {
		return controller.getUid();
	}
	
	
	/**
	 * Getter for the spectrumIdentification object of the PRIDE implementation
	 * 
	 * @return
	 */
	private SpectrumIdentification getSpectrumIdentification() {
		return controller.getPeptideByIndex(proteinID, peptideID).getSpectrumIdentification();
	}
	
	
	@Override
	public Double getScore(String scoreAccession) {
		if ((cachedScoreAccession != null) && (cachedScoreAccession.equals(scoreAccession))) {
			return cachedScore;
		} else if (CvScore.PSI_PSM_LEVEL_FDRSCORE.getAccession().equals(scoreAccession)) {
			return getFDRScore();
		} else if (CvScore.PSI_PSM_LEVEL_LOCAL_FDR.getAccession().equals(scoreAccession)) {
			return getFDR();
		} else if (CvScore.PSI_PSM_LEVEL_Q_VALUE.getAccession().equals(scoreAccession)) {
			return getQValue();
		} else {
			CvTermReference cvTermRef = CvTermReference.getCvRefByAccession(scoreAccession);
			if (cvTermRef != null) {
				List<Number> scores = 
						getSpectrumIdentification().getScore().getScores(cvTermRef);
				
				if (scores.size() > 0) {
					return scores.get(0).doubleValue();
				}
			}
			
			return null;
		}
	}
	
	
	@Override
	public List<String> getBaseScoreAccessions() {
		List<String> scoreAccessions = new ArrayList<String>();
		
		for (CvTermReference cvTerm
				: getSpectrumIdentification().getScore().getCvTermReferenceWithValues()) {
			scoreAccessions.add(cvTerm.getAccession());
		}
		
		return scoreAccessions;
	}
	
	
	@Override
	public String getSpectrumId() {
		// not supported by pride
		return null;
	}
	
	
	@Override
	public Double getExperimentalMassToCharge() {
		return getSpectrumIdentification().getExperimentalMassToCharge();
	}
	
	
	@Override
	public Double getDeltaMass() {
		return getSpectrumIdentification().getExperimentalMassToCharge() - 
				getSpectrumIdentification().getCalculatedMassToCharge();
	}
	
	
	@Override
	public Double getRetentionTime() {
		// not supported by pride
		return null;
	}


	@Override
	public Integer getCharge() {
		return getSpectrumIdentification().getChargeState();
	}
	
	
	@Override
	public Integer getMissedCleavages() {
		// missed cleavages are not given in mzIdentML etc.
		return null;
	}
	
	
	@Override
	public String getSequence() {
		return controller.getPeptideByIndex(proteinID, peptideID).getSequence();
	}
	
	
	@Override
	public List<Modification> getModifications() {
		List<Modification> modifications = new ArrayList<Modification>();
		
		for (uk.ac.ebi.pride.utilities.data.core.Modification prideMod
				: controller.getPeptideByIndex(proteinID, peptideID).getModifications()) {
			Modification mod = new Modification();
			
			Double massDelta = null;
			int count = 0;
			for (Double delta :	prideMod.getAvgMassDelta()) {
				if (delta != null) {
					if (massDelta != null) {
						massDelta += delta;
					} else {
						massDelta = delta;
					}
					count++;
				}
			}
			if (massDelta != null) {
				massDelta /= count;
				mod.setAvgMassDelta(massDelta);
			}
			
			mod.setLocation(prideMod.getLocation());
			
			massDelta = null;
			count = 0;
			for (Double delta :	prideMod.getMonoisotopicMassDelta()) {
				if (delta != null) {
					if (massDelta != null) {
						massDelta += delta;
					} else {
						massDelta = delta;
					}
					count++;
				}
			}
			if (massDelta != null) {
				massDelta /= count;
				mod.setMonoisotopicMassDelta(massDelta);
			}
			
			mod.getResidues().addAll(prideMod.getResidues());
			
			mod.getCvParam().addAll(PrideUtilities.convertCvParams(prideMod.getCvParams()));
			
			modifications.add(mod);
		}
		
		return modifications;
	}
	
	
	@Override
	public List<AbstractParam> getParams() {
		List<AbstractParam> params = new ArrayList<AbstractParam>();
		
		params.addAll(PrideUtilities.convertCvParams(getSpectrumIdentification().getCvParams()));
		params.addAll(PrideUtilities.convertUserParams(getSpectrumIdentification().getUserParams()));
		
		return params;
	}
}
