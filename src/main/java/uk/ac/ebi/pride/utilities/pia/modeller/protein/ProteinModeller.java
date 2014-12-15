package uk.ac.ebi.pride.utilities.pia.modeller.protein;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateStructure;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.protein.inference.AbstractProteinInference;
import uk.ac.ebi.pride.utilities.pia.modeller.protein.inference.InferenceProteinGroup;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide.PeptideScoring;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.protein.ProteinScoring;


/**
 * This class handles the proteins and their inference in a given
 * {@link uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateStructure}.
 *
 * @author julian
 *
 */
public class ProteinModeller {

    /** logger for this class */
    private static final Logger logger =  LoggerFactory.getLogger(ProteinModeller.class);

    /** the associated intermediate structure */
    private IntermediateStructure intermediateStructure;

    /** the maximal number of allowed threads */
    private int allowedThreads;

    /** the inferred proteins */
    private List<InferenceProteinGroup> inferredProteins;

    /** the peptide scoring used for protein inference */
    private PeptideScoring usedPeptideScoring;

    /** the protein scoring used for protein inference */
    private ProteinScoring usedProteinScoring;

    /** the used filters */
    private List<AbstractFilter> usedFilters;

    /** whether modifications were considered */
    private Boolean usedConsiderModifications;




    public ProteinModeller(IntermediateStructure intermediateStructure, int allowedThreads) {
        this.intermediateStructure = intermediateStructure;
        this.allowedThreads = allowedThreads;
        this.inferredProteins = null;

        this.usedPeptideScoring = null;
        this.usedProteinScoring = null;
        this.usedFilters = null;
        this.usedConsiderModifications = null;
    }


    /**
     * Runs the protein inference with the given settings and stores the protein
     * groups in a list, which can be returned by {@link #getInferredProteins()}
     *
     * @param peptideScoring
     * @param proteinScoring
     * @param proteinInferenceClass
     * @param filters
     * @param considerModifications
     */
    public void infereProteins(PeptideScoring peptideScoring, ProteinScoring proteinScoring,
                               Class<? extends AbstractProteinInference> proteinInferenceClass, List<AbstractFilter> filters,
                               boolean considerModifications) {
        AbstractProteinInference proteinInference;

        try {
            Constructor<? extends AbstractProteinInference> constructor =
                    proteinInferenceClass.getConstructor(
                            IntermediateStructure.class, PeptideScoring.class, ProteinScoring.class, List.class, Integer.class);

            proteinInference = constructor.newInstance(intermediateStructure, peptideScoring, proteinScoring, filters, allowedThreads);
        } catch (Exception ex) {
            logger.error("Could not initialize protein inference for " + proteinInferenceClass.getCanonicalName(), ex);
            inferredProteins = null;
            return;
        }

        inferredProteins = proteinInference.calculateInference(considerModifications);

        usedPeptideScoring = peptideScoring;
        usedProteinScoring = proteinScoring;
        if (filters != null) {
            usedFilters = new ArrayList<AbstractFilter>(filters);
        }
        usedConsiderModifications = considerModifications;
    }


    /**
     * Returns a list of the inferred protein groups. These protein groups must
     * be calculated by {@link #infereProteins(PeptideScoring, ProteinScoring, Class, List, boolean)}
     * beforehand.
     *
     * @return
     */
    public List<InferenceProteinGroup> getInferredProteins() {
        if (inferredProteins != null) {
            return new ArrayList<InferenceProteinGroup>(inferredProteins);
        } else {
            return null;
        }
    }
}
