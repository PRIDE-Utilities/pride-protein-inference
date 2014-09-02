package uk.ac.ebi.pride.utilities.pia.modeller.protein.inference;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import uk.ac.ebi.pride.utilities.data.core.DBSequence;
import uk.ac.ebi.pride.utilities.data.core.Peptide;
import uk.ac.ebi.pride.utilities.data.core.PeptideEvidence;
import uk.ac.ebi.pride.utilities.data.core.Protein;
import uk.ac.ebi.pride.utilities.data.core.ProteinGroup;
import uk.ac.ebi.pride.utilities.data.core.Score;
import uk.ac.ebi.pride.utilities.data.core.SpectrumIdentification;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptide;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateProtein;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.ScoringItemType;

/**
 * This group bundles the information about an inferred protein group.
 * 
 * @author julian
 *
 */
public class InferenceProteinGroup {
	
	/** a unique identifier */
	private String ID;
	
	/** the proteins in this group */
	private Set<IntermediateProtein> proteins;
	
	/** the identified peptides leading to this group */
	private Map<Comparable, IntermediatePeptide> intermediatePeptides;
	
	/** any subgroups */
	private Set<InferenceProteinGroup> subGroups;
	
	/** whether modifications were considered to distinguish pepttides of this protein group */
	boolean considerModifications;
	
	/** the protein score */
	private Double score;
	
	/** The IDs of the intermediatePeptides mapping to the type of scoring. If a
	 *  peptide's ID is not a key in the map, it is assumed to not score */
	private Map<Comparable, ScoringItemType> peptideScorings;
	
	
	/**
	 * Basic constructor
	 *
	 */
	public InferenceProteinGroup(String id, boolean considerModifications) {
		this.ID = id;
		this.proteins = new HashSet<IntermediateProtein>();
		this.intermediatePeptides = new HashMap<Comparable, IntermediatePeptide>();
		this.subGroups = new HashSet<InferenceProteinGroup>();
		this.considerModifications = considerModifications;
		this.score = Double.NaN;
		this.peptideScorings = null;
	}
	
	
	public String getID() {
		return ID;
	}
	
	
	/**
	 * Adds one protein to the set of proteins
	 * @param protein
	 * @return
	 */
	public boolean addProtein(IntermediateProtein protein) {
		return proteins.add(protein);
	}
	
	
	/**
	 * Returns the set of proteins
	 * @return
	 */
	public Set<IntermediateProtein> getProteins() {
		return proteins;
	}
	
	
	/**
	 * Adds one IntermediatePeptide to the map of peptides
	 * @return
	 */
	public IntermediatePeptide addPeptide(IntermediatePeptide peptide) {
		Comparable peptideID =
				AbstractProteinInference.getPeptideKey(peptide, considerModifications);
		return intermediatePeptides.put(peptideID, peptide);
	}
	
	
	/**
	 * Adds a collection of IntermediatePeptides to the map of peptides
	 * @return
	 */
	public boolean addPeptides(Collection<IntermediatePeptide> peptides) {
		boolean changes = false;
		
		for (IntermediatePeptide peptide : peptides) {
			Comparable peptideID =
					AbstractProteinInference.getPeptideKey(peptide, considerModifications);
			changes |= (intermediatePeptides.put(peptideID, peptide) != null);
		}
		
		return changes;
	}
	
	
	/**
	 * Returns the set of IntermediatePeptides
	 * @return
	 */
	public Set<IntermediatePeptide> getPeptides() {
		return new HashSet<IntermediatePeptide>(intermediatePeptides.values());
	}
	
	
	/**
	 * Getter for the protein score. If the score is not given, it may be null
	 * 
	 * @return
	 */
	public Double getScore() {
		return score;
	}
	
	
	/**
	 * Sets the score of the protein.
	 * 
	 * @param score
	 */
	public void setScore(Double score) {
		this.score = score;
	}
	
	
	/**
	 * Adds a protein group to the subGroups of this group
	 * @param proteinGroup
	 * @return
	 */
	public boolean addSubgroup(InferenceProteinGroup proteinGroup) {
		return subGroups.add(proteinGroup);
	}
	
	
	/**
	 * Returns the set of sub-groups
	 * @return
	 */
	public Set<InferenceProteinGroup> getSubGroups() {
		return subGroups;
	}
	
	
	
	/**
	 * Sets the scoring type of the stated peptide to the given {@link ScoringItemType}.
	 *
	 */
	public void setPeptidesScoringType(IntermediatePeptide peptide, ScoringItemType type) {
		if (peptideScorings == null) {
			peptideScorings = new HashMap<Comparable, ScoringItemType>();
		}
		
		Comparable peptideID =
				AbstractProteinInference.getPeptideKey(peptide, considerModifications);
		if (intermediatePeptides.containsKey(peptideID)) {
			peptideScorings.put(peptideID, type);
		}
	}
	
	
	/**
	 * Returns the {@link ScoringItemType} of the peptide
	 * @return
	 */
	public ScoringItemType getPeptidesScoringType(IntermediatePeptide peptide) {
		Comparable peptideID =
				AbstractProteinInference.getPeptideKey(peptide, considerModifications);
		
		if ((peptideScorings == null) || !peptideScorings.containsKey(peptideID)) {
			return ScoringItemType.NOT_SCORING;
		} else {
			return peptideScorings.get(peptideID);
		}
	}
	
	
	/**
	 * Removes all information about which peptides were used for scoring.
	 */
	public void removeAllScoringInformation() {
		peptideScorings = null;
	}
	
	
	/**
	 * Returns the peptides with the given {@link ScoringItemType}, which also pass the 
	 * @return
	 */
	public Set<IntermediatePeptide> getPeptidesWithScoringType(ScoringItemType type) {
		HashSet<IntermediatePeptide> peptides =  new HashSet<IntermediatePeptide>();
		
		if ((peptideScorings == null) && (ScoringItemType.NOT_SCORING.equals(type))) {
			return getPeptides();
		} else {
			for (IntermediatePeptide peptide : getPeptides()) {
				if (getPeptidesScoringType(peptide).equals(type)) {
					peptides.add(peptide);
				}
			}
		}
		
		return peptides;
	}
	
	
	/**
	 * Creates a {@link ProteinGroup} (a protein ambiguity group in mzIdentML)
	 * from this inferred protein group.
	 * @return
	 */
	public ProteinGroup createProteinGroup() {
		Set<DBSequence> dbSequences = new HashSet<DBSequence>();
		for (IntermediateProtein interProt : proteins) {
			dbSequences.add(interProt.getDBSequence());
		}
		
		Set<Peptide> peptideHypotheses = new HashSet<Peptide>();
		for (IntermediatePeptide peptide : getPeptides()) {
			for (IntermediatePeptideSpectrumMatch psm : peptide.getPeptideSpectrumMatches()) {
				SpectrumIdentification specID = psm.getSpectrumIdentification();
				for (PeptideEvidence pepEvidence : specID.getPeptideEvidenceList()) {
					Peptide pep = new Peptide(pepEvidence, specID);
					peptideHypotheses.add(pep);
					
					// TODO: set the scoring information
				}
			}
		}
		
		List<Protein> proteins = new ArrayList<Protein>();
		for (DBSequence dbSeq : dbSequences) {
			Score score = null;
			double sequenceCoverage = -1; // TODO: include calculation of this coverage
			
			Protein protein = new Protein(
					ID + "_" + dbSeq.getAccession(),
					dbSeq.getAccession(),
					dbSeq,
					true /*passThreshold*/,
					new ArrayList<Peptide>(peptideHypotheses),
					score,
					-1 /*threshold*/,
					sequenceCoverage,
					null /*gel*/);
			
			
			proteins.add(protein);
		}
		
		// the sub-groups proteins will become proteins which do not pass the threshold
		for (InferenceProteinGroup subGroup : subGroups) {
			ProteinGroup subProteinGroup = subGroup.createProteinGroup();
			
			for (Protein subProtein : subProteinGroup.getProteinDetectionHypothesis()) {
				subProtein.setId(ID + "_" + subProtein.getDbSequence().getAccession());
				subProtein.setPassThreshold(false);
				proteins.add(subProtein);
			}
		}
		
		return new ProteinGroup(ID, ID, proteins);
	}
	
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
		
		InferenceProteinGroup group = (InferenceProteinGroup)obj;

        return ID == group.ID && !((proteins != null) ? !proteins.equals(group.proteins) : (group.proteins != null)) && !((intermediatePeptides != null) ? !intermediatePeptides.equals(group.intermediatePeptides) : (group.intermediatePeptides != null)) && !((subGroups != null) ? !subGroups.equals(group.subGroups) : (group.subGroups != null));
    }
	
	
	@Override
	public int hashCode() {
        int result = ID.hashCode();
        result = 31 * result + ((proteins != null) ? proteins.hashCode() : 0);
        result = 31 * result + ((intermediatePeptides != null) ? intermediatePeptides.hashCode() : 0);
        result = 31 * result + ((subGroups != null) ? subGroups.hashCode() : 0);
        return result;
	}
}
