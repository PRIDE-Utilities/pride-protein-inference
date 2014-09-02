package uk.ac.ebi.pride.utilities.pia.modeller.filter.protein;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptide;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterComparator;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterType;
import uk.ac.ebi.pride.utilities.pia.modeller.protein.inference.InferenceProteinGroup;

/**
 * Filters for the number of PSMs on protein level.
 * <p>
 * If filters are set on the peptides, only the passing PSMs are counted.
 * 
 * @author julian
 *
 */
public class ProteinNrPSMsFilter extends AbstractFilter {
	
	protected static final String shortName = "protein_nr_psms_filter";
	
	private static final String name = "Number of PSMs Filter for Protein";
	
	private static final String filteringName = "#PSMs (Protein)";
	
	private static final FilterType filterType = FilterType.numerical;
	
	private Integer value;
	
	
	public ProteinNrPSMsFilter(FilterComparator arg, Integer value, boolean negate) {
		this.comparator = arg;
		this.value = value;
		this.negate = negate;
	}
	
	
	@Override
	public String getShortName() {
		return shortName;
	}
	
	@Override
	public String getName() {
		return name;
	}
	
	@Override
	public String getFilteringName() {
		return filteringName;
	}
	
	@Override
	public Object getFilterValue() {
		return value;
	}
	
	@Override
	public FilterType getFilterType() {
		return filterType;
	}
	
	@Override
	public Object getObjectsValue(Object o) {
		if (o instanceof InferenceProteinGroup) {
			int nrPSMs = 0;
			
			for (IntermediatePeptide peptide : ((InferenceProteinGroup) o).getPeptides()) {
				nrPSMs += peptide.getNumberOfPeptideSpectrumMatches();
			}
			
			return nrPSMs;
		} else {
			return null;
		}
	}
	
	@Override
	public boolean supportsClass(Object c) {
        return c instanceof InferenceProteinGroup;
	}
}
