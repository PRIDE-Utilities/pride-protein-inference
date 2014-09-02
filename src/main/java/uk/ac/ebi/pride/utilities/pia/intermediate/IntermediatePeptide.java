package uk.ac.ebi.pride.utilities.pia.intermediate;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterUtilities;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.ScoringItemType;


/**
 * A peptide class, which is used in the intermediate structure.
 * 
 * @author julian
 * @author ypriverol
 *
 */
public class IntermediatePeptide {
	
	/** Sequence of the peptide */
	private String sequence;
	
	/** the connected group of this peptides */
	private IntermediateGroup group;
	
	/** List of identifications for this peptide */
	private Map<Comparable, IntermediatePeptideSpectrumMatch> peptideSpectrumMatches;
	
	/** The PSMs passing a given filter */
	private Set<Comparable> psmsPassingFilter;
	
	/** the peptide score */
	private Double score;
	
	/** The IDs of the PSMs mapping to the type of scoring. If a PSM's ID is not
	 *  a key in the map, it is assumed to not score */
	private Map<Comparable, ScoringItemType> psmScorings;
	
	
	/**
	 * Basic constructor, only initializes the sequence
	 *  
	 * @param sequence
	 */
	public IntermediatePeptide(String sequence) {
		this.sequence = sequence;
		this.group = null;
		this.peptideSpectrumMatches = new HashMap<Comparable, IntermediatePeptideSpectrumMatch>();
		this.psmsPassingFilter = null;
		this.score = Double.NaN;
		this.psmScorings = null;
	}
	
	
	/**
	 * Returns an identifier for the peptide
	 * 
	 * @return
	 */
	public Comparable getID() {
		return computeID(sequence);
	}
	
	
	/**
	 * Returns the ID comoputed by the given arguments.
	 * <p>
	 * This is for probable future integration of modifications
	 * 
	 * @param sequence
	 * @return
	 */
	public static String computeID(String sequence) {
		return sequence;
	}
	
	
	/**
	 * Getter for the sequence.
	 * 
	 * @return
	 */
	public String getSequence() {
		return sequence;
	}
	
	
	/**
	 * sets the group of this peptide
	 * @param group
	 */
	public void setGroup(IntermediateGroup group) {
		this.group = group;
	}
	
	
	/**
	 * returns the group of this peptide
	 */
	public IntermediateGroup getGroup() {
		return group;
	}
	
	
	/**
	 * Adds the given spectrumIdentification to the peptide's PSMs
	 * <p>
	 * Is the PSM is added, it is also added to the filter passing PSMs.
	 * 
	 * @param spectrum
	 * @return true if this peptide did not already contain the spectrum identification
	 */
	public boolean addPeptideSpectrumMatch(IntermediatePeptideSpectrumMatch spectrumIdentification) {
		if (!peptideSpectrumMatches.containsKey(spectrumIdentification.getID())) {
			peptideSpectrumMatches.put(spectrumIdentification.getID(), spectrumIdentification);
			
			if (psmsPassingFilter != null) {
				psmsPassingFilter.add(spectrumIdentification.getID());
			}
			return true;
		}
		
		return false;
	}
	
	
	/**
	 * Getter for the PSMs.
	 * <p>
	 * If a passing set is created (i.e. data is filtered), only the passing
	 * PSMs are returned. 
	 * 
	 * @return
	 */
	public List<IntermediatePeptideSpectrumMatch> getPeptideSpectrumMatches() {
		if (psmsPassingFilter != null) {
			List<IntermediatePeptideSpectrumMatch> psms =
					new ArrayList<IntermediatePeptideSpectrumMatch>(psmsPassingFilter.size());
			for (Map.Entry<Comparable, IntermediatePeptideSpectrumMatch> psmIt : peptideSpectrumMatches.entrySet()) {
				if (psmsPassingFilter.contains(psmIt.getKey())) {
					psms.add(psmIt.getValue());
				}
			}
			return psms;
		} else {
			return new ArrayList<IntermediatePeptideSpectrumMatch>(peptideSpectrumMatches.values());
		}
	}
	
	
	/**
	 * Getter for the number of PSMs.
	 * <p>
	 * If a passing set is created (i.e. data is filtered), only the passing
	 * PSMs are returned.
	 *  
	 * @return
	 */
	public Integer getNumberOfPeptideSpectrumMatches() {
		if (psmsPassingFilter != null) {
			return psmsPassingFilter.size();
		} else {
			return peptideSpectrumMatches.size();
		}
	}
	
	
	/**
	 * Getter for the PSMs, also the ones, which do not pass the filters.
	 * @return
	 */
	public List<IntermediatePeptideSpectrumMatch> getAllPeptideSpectrumMatches() {
		return new ArrayList<IntermediatePeptideSpectrumMatch>(peptideSpectrumMatches.values());
	}
	
	
	/**
	 * Filters the PSMs using the given filters. Any prior filtering will be
	 * deleted.
	 * 
	 * @param filters
	 */
	public void filterPSMs(List<AbstractFilter> filters) {
		psmsPassingFilter = new HashSet<Comparable>();
		
		for (IntermediatePeptideSpectrumMatch psm : getAllPeptideSpectrumMatches()) {
			if (FilterUtilities.satisfiesFilterList(psm, filters)) {
				psmsPassingFilter.add(psm.getID());
				
			}
		}
	}
	
	
	/**
	 * Getter for the peptide score. If the score is not given, it may be null
	 * or {@value Double#NaN}. 
	 * 
	 * @return
	 */
	public Double getScore() {
		return score;
	}
	
	
	/**
	 * Sets the score of the peptide.
	 * 
	 * @param score
	 */
	public void setScore(Double score) {
		this.score = score;
	}
	
	
	/**
	 * Sets the scoring type of the stated PSM to the given {@link ScoringItemType}.
	 * 
	 * @param psmID
	 */
	public void setPSMsScoringType(IntermediatePeptideSpectrumMatch psm, ScoringItemType type) {
		if (psmScorings == null) {
			psmScorings = new HashMap<Comparable, ScoringItemType>();
		}
		
		Comparable psmID = psm.getID();
		if (peptideSpectrumMatches.containsKey(psmID)) {
			psmScorings.put(psmID, type);
		}
	}
	
	
	/**
	 * Returns the {@link ScoringItemType} of the PSM
	 * 
	 * @param psmID
	 * @return
	 */
	public ScoringItemType getPSMsScoringType(IntermediatePeptideSpectrumMatch psm) {
		Comparable psmID = psm.getID();
		
		if ((psmScorings == null) || !psmScorings.containsKey(psmID)) {
			return ScoringItemType.NOT_SCORING;
		} else {
			return psmScorings.get(psmID);
		}
	}
	
	
	/**
	 * Removes all information about which PSMs were used for scoring.
	 */
	public void removeAllScoringInformation() {
		psmScorings = null;
	}
	
	
	/**
	 * Returns the PSMs with the given {@link ScoringItemType}, which also pass
	 * any filterings.
	 * 
	 * @return
	 */
	public Set<IntermediatePeptideSpectrumMatch> getPSMsWithScoringType(ScoringItemType type) {
		HashSet<IntermediatePeptideSpectrumMatch> psms =  new HashSet<IntermediatePeptideSpectrumMatch>();
		
		if ((psmScorings == null) && (ScoringItemType.NOT_SCORING.equals(type))) {
			return new HashSet<IntermediatePeptideSpectrumMatch>(getPeptideSpectrumMatches());
		} else {
			for (IntermediatePeptideSpectrumMatch psm : getPeptideSpectrumMatches()) {
				if (getPSMsScoringType(psm).equals(type)) {
					psms.add(psm);
				}
			}
		}
		
		return psms;
	}
	
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
        if (obj == null || !(obj instanceof IntermediatePeptide)) return false;
		
		IntermediatePeptide peptide = (IntermediatePeptide)obj;
		
		if (!sequence.equals(peptide.getSequence())) return false;
		if (psmsPassingFilter != null) {
			if (!psmsPassingFilter.equals(peptide.psmsPassingFilter)) {
				return false;
			}
		} else {
			if (peptide.psmsPassingFilter != null) {
				return false;
			}
		}
		
		return !((group != null) ? !group.getID().equals(peptide.getGroup().getID()) : (peptide.getGroup() != null)) &&
				!((peptideSpectrumMatches != null) ? !peptideSpectrumMatches.equals(peptide.peptideSpectrumMatches) : (peptide.peptideSpectrumMatches != null));
	}
	
	
	@Override
	public int hashCode() {
        int result = (sequence != null) ? sequence.hashCode() : 0;
        result = 31 * result + ((group != null) ? group.getID().hashCode() : 0);
        result = 31 * result + ((peptideSpectrumMatches != null) ? peptideSpectrumMatches.hashCode() : 0);
        return result;
	}
}
