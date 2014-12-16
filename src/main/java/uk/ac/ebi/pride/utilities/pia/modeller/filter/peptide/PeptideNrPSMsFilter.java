package uk.ac.ebi.pride.utilities.pia.modeller.filter.peptide;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptide;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterComparator;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterType;


/**
 * Filters for the number of PSMs on peptide level.
 * <p>
 * If filters are set on the peptide, only the passing PSMs are counted.
 * 
 * @author julian
 *
 */
public class PeptideNrPSMsFilter extends AbstractFilter {
	
	protected static final String shortName = "peptide_nr_psms_filter";
	
	private static final String name = "Number of PSMs Filter for Peptide";
	
	private static final String filteringName = "#PSMs (Peptide)";
	
	private static final FilterType filterType = FilterType.numerical;
	
	private Integer value;
	
	
	public PeptideNrPSMsFilter(FilterComparator arg, Integer value, boolean negate) {
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
		if (o instanceof IntermediatePeptide) {
			return ((IntermediatePeptide) o).getNumberOfPeptideSpectrumMatches();
		} else {
			return null;
		}
	}

	@Override
	public boolean supportsClass(Object c) {
		if (c instanceof IntermediatePeptide) {
			return true;
		} else {
			return false;
		}
	}

}
