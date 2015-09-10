package uk.ac.ebi.pride.utilities.pia.modeller.filter.protein;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateProtein;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterComparator;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterType;


/**
 * Filters for accessions on PSM level.
 * 
 * @author julian
 *
 */
public class ProteinAccessionFilter extends AbstractFilter {
	
	protected static final String shortName = "protein_accession_filter";
	
	private static final String name = "Accession Filter for a single protein";
	
	private static final String filteringName = "Accession (protein)";
	
	private static final FilterType filterType = FilterType.literal;
	
	private String value;
	
	
	public ProteinAccessionFilter(FilterComparator arg, String value, boolean negate) {
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
		if (o instanceof IntermediateProtein) {
			return ((IntermediateProtein) o).getAccession();
		} else {
			return null;
		}
	}

	@Override
	public boolean supportsClass(Object c) {
        return c instanceof IntermediateProtein;
	}

}
