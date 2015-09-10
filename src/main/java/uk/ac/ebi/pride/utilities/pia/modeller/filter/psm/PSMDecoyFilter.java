package uk.ac.ebi.pride.utilities.pia.modeller.filter.psm;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterComparator;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterType;


/**
 * Filter for decoys on PSM level.
 * 
 * @author julian
 *
 */
public class PSMDecoyFilter extends AbstractFilter {
	
	protected static final String shortName = "psm_decoy_filter";
	
	private static final String name = "Decoy Filter for PSM";
	
	private static final String filteringName = "Decoy (PSM)";
	
	private static final FilterType filterType = FilterType.bool;
	
	private Boolean value;
	
	
	public PSMDecoyFilter(FilterComparator arg, Boolean value, boolean negate) {
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
		if (o instanceof IntermediatePeptideSpectrumMatch) {
			return ((IntermediatePeptideSpectrumMatch) o).getIsDecoy();
		} else {
			return null;
		}
	}
	
	
	@Override
	public boolean supportsClass(Object c) {
        return (c instanceof IntermediatePeptideSpectrumMatch);
	}

}
