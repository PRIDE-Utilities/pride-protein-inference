package uk.ac.ebi.pride.utilities.pia.modeller.fdr;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateProtein;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.protein.ProteinAccessionFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.ScoreUtilities;


/**
 * This class provides functionality for the FDR estimation using decoy proteins
 *
 * @author julian
 *
 */
public class FDRUtilities {

    /** logger for this class */
    private static final Logger logger =  LoggerFactory.getLogger(FDRUtilities.class);


    /**
     * Calculate the FDR on the given score sorted List of
     * {@link FDRComputableByDecoys} objects.
     *
     * The item with the best score must be the first in the list, the worst
     * score the last.
     *
     */
    public static <T extends FDRComputableByDecoys> void calculateFDR(List<T> items, String scoreAccession) {
        double fdr;

        Double rankScore = Double.NaN;
        List<T> rankItems = new ArrayList<T>();

        int nrTargets = 0;
        int nrDecoys = 0;

        for (T item : items) {
            Double itemScore = item.getScore(scoreAccession);
            if (!rankScore.equals(itemScore)) {
                // this is a new rank, calculate FDR
                if (!rankScore.equals(Double.NaN) && (nrTargets < 1)) {
                    // only decoys until now -> set FDR to infinity
                    fdr = Double.POSITIVE_INFINITY;
                } else {
                    fdr = (double)nrDecoys / nrTargets;
                }

                for (T rankItem : rankItems) {
                    rankItem.setFDR(fdr);
                }

                rankScore = itemScore;
                rankItems.clear();
            }

            // check for decoy
            if (item.getIsDecoy()) {
                nrDecoys++;
            } else {
                nrTargets++;
            }

            rankItems.add(item);
        }

        // calculate the last rank
        if (nrTargets < 1) {
            // only decoys until now -> set FDR to infinity
            fdr = Double.POSITIVE_INFINITY;
        } else {
            fdr = (double)nrDecoys / nrTargets;
        }

        for (T rankItem : rankItems) {
            rankItem.setFDR(fdr);
        }


        // at last calculate the q-values
        // for this, iterate backwards through the list
        ListIterator<T> it  = items.listIterator(items.size());
        Double qValue = Double.NaN;

        while (it.hasPrevious()) {
            T item = it.previous();

            if ((qValue.compareTo(Double.NaN) == 0) ||
                    (item.getFDR() < qValue)) {
                qValue = item.getFDR();
            }

            item.setQValue(qValue);
        }
    }


    /**
     * Calculates the FDR score of the report. To do this, the report must have
     * FDR values and be sorted.
     *
     * @param items the list of items, for which the FDR will be calculated
     * @param scoreAccession the accession of the score used for FDR calculation
     * @param oboLookup whether ob olookup should be performed, if the score is not hard-coded
     */
    public static <T extends FDRComputableByDecoys> void calculateFDRScore(
            List<T> items, String scoreAccession, boolean oboLookup) {
        if (items.size() < 2) {
            // no calculation for empty list possible
            return;
        }

        // get the stepPoints of the q-values
        List<Integer> stepPoints = new ArrayList<Integer>();
        ListIterator<T> it = items.listIterator(items.size());
        double qValue = Double.NaN;
        int nrDecoys = 0;
        int nrTargets = 0;

        while (it.hasPrevious()) {
            T item = it.previous();

            if ((Double.compare(qValue, Double.NaN) != 0) &&
                    (item.getQValue() < qValue)) {
                stepPoints.add(it.nextIndex()+1);
            }

            qValue = item.getQValue();

            if (item.getIsDecoy()) {
                nrDecoys++;
            } else {
                nrTargets++;
            }
        }

        // calculate the FDR scores
        double g;
        double qLast, qNext;
        double sLast, sNext;

        Collections.sort(stepPoints);
        ListIterator<Integer> stepIterator = stepPoints.listIterator();
        Integer nextStep;

        if (ScoreUtilities.isHigherScoreBetter(scoreAccession, oboLookup)) {
            // get the score of the first entry + (difference between first entry and first decoy) / (index of first decoy)  (to avoid FDRScore = 0)
            Double zeroScore = items.get(0).getScore(scoreAccession);
            if (stepPoints.size() > 0) {
                sLast = zeroScore +
                        (zeroScore - items.get(stepPoints.get(0)).getScore(scoreAccession)) / stepPoints.get(0);
            } else {
                sLast = zeroScore +
                        (zeroScore - items.get(items.size()-1).getScore(scoreAccession)) / items.size()-1;
            }
        } else {
            // or 0, if not higherscorebetter
            sLast = 0;
        }
        qLast = 0;

        if (stepIterator.hasNext()) {
            nextStep = stepIterator.next();

            sNext = items.get(nextStep).getScore(scoreAccession);
            qNext = items.get(nextStep).getQValue();
        } else {
            // we add an artificial decoy to the end...
            nextStep = items.size();

            sNext = items.get(items.size()-1).getScore(scoreAccession);
            qNext = (nrTargets == 0) ? Double.POSITIVE_INFINITY : (double)(nrDecoys + 1) / nrTargets;
        }

        g = (qNext-qLast) / (sNext-sLast);

        it = items.listIterator();
        while (it.hasNext()) {
            T item = it.next();

            if (nextStep == it.nextIndex()-1) {
                if (stepIterator.hasNext()) {
                    sLast = sNext;
                    qLast = qNext;

                    nextStep = stepIterator.next();

                    sNext = items.get(nextStep).getScore(scoreAccession);
                    qNext = items.get(nextStep).getQValue();
                }

                g = (qNext-qLast) / (sNext-sLast);
            }

            item.setFDRScore((item.getScore(scoreAccession)-sLast)*g + qLast);
        }
    }


    /**
     * This function marks the intermediateProteins of the PSMs, which pass the
     * given decoysFilter, as decoys for a subsequent FDR estimation.
     *
     * @param psms the PSMs, of which the  proteins are tagged
     * @return the number of decoys in the list
     */
    public static int markDecoys(List<IntermediatePeptideSpectrumMatch> psms,
                                 ProteinAccessionFilter decoysFilter) {
        int count = 0;
        Set<String> proteinsDone = new HashSet<String>();

        for (IntermediatePeptideSpectrumMatch psm : psms) {
            for (IntermediateProtein protein : psm.getPeptide().getAllProteins()) {
                String protAccession = protein.getAccession();
                if (!proteinsDone.contains(protAccession)) {
                    if (decoysFilter.satisfiesFilter(protein)) {
                        protein.setIsDecoy(true);
                        count++;
                    } else {
                        protein.setIsDecoy(false);
                    }
                    proteinsDone.add(protAccession);
                }
            }
        }

        logger.info("decoys marked (" + count + "/" + proteinsDone.size() + ")");
        return count;
    }
}