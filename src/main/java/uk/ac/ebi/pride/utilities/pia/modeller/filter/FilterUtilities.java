package uk.ac.ebi.pride.utilities.pia.modeller.filter;

import java.util.ArrayList;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.pride.utilities.pia.modeller.filter.psm.PSMScoreFilter;


/**
 * Some helper methods for the filtering.
 * 
 * @author julian
 *
 */
public final class FilterUtilities {
	public static final Logger logger = LoggerFactory.getLogger(FilterUtilities.class);

    /**
     * Do not instantiate this class, all methods are static!
     */
    private FilterUtilities() {
    	throw new AssertionError();
    }
    
    
	/**
	 * Checks whether all the filters in the given List are satisfied for the
	 * given object.
	 * 
	 * @param item
	 * @param fileID
	 * @param filters
	 * @return
	 */
	public static boolean satisfiesFilterList(Object item,
			List<AbstractFilter> filters) {
		boolean satisfiesAllFilters = true;
		
		for (AbstractFilter filter : filters) {
			if (filter.supportsClass(item)) {
				satisfiesAllFilters &= filter.satisfiesFilter(item);
			}
		}
		
		return satisfiesAllFilters;
	}
    
	
    /**
	 * gets the available comparators for the filter given by the short name
	 * 
	 * @param filterShort
	 * @return
	 */
	public static List<FilterComparator> getAvailableComparators(String filterShort) {
		if (filterShort != null) {
			FilterType type = getFilterType(filterShort);
			if (type != null) {
				return type.getAvailableComparators();
			}
		}
		
		return new ArrayList<FilterComparator>();
	}
	
	
	/**
	 * Returns the FilterType of the given filter.
	 * 
	 * @param filterShort
	 * @return
	 */
	public static FilterType getFilterType(String filterShort) {
		
		if (filterShort != null) {
			if (filterShort.startsWith(PSMScoreFilter.prefix)) {
				return PSMScoreFilter.filterType;
			}/* else if (filterShort.startsWith(PeptideScoreFilter.prefix)) {
				return PeptideScoreFilter.filterType;
			} else if (filterShort.equals(ProteinScoreFilter.shortName)) {
				return ProteinScoreFilter.filterType;
			} else if (filterShort.startsWith(PSMTopIdentificationFilter.prefix)) {
				return PSMTopIdentificationFilter.filterType;
			} else {
				
				for (RegisteredFilters filter : RegisteredFilters.values()) {
					if (filter.getShortName().equals(filterShort)) {
						return filter.getFilterType();
					}
				}
			}*/
		}
		
		return null;
	}
}
