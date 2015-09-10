package uk.ac.ebi.pride.utilities.pia.modeller.filter.psm;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterComparator;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterType;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.CvScore;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.ScoreUtilities;



/**
 * Filters for a score on PSM level.
 * 
 * @author julian
 *
 */
public class PSMScoreFilter extends AbstractFilter {
	
	protected final String shortName;
	
	protected final String name;
	
	protected final String filteringName;
	
	public static final FilterType filterType = FilterType.numerical;
	
	private Double value;
	
	public static String prefix = "psm_score_filter_";
	
	private String cvAccession;
	
	
	public PSMScoreFilter(FilterComparator arg, Double value, boolean negate,
			String scoreAccession, boolean oboLookup) {
		cvAccession = null;
		
		CvScore cvScore = CvScore.getCvRefByAccession(scoreAccession);
		if (cvScore == null) {
			// try with "MS" enhenced name
			cvScore = CvScore.getCvRefByAccession("MS:" + scoreAccession);
		}
		if (cvScore != null) {
			cvAccession = cvScore.getAccession();
		} else if (oboLookup) {
			cvAccession = ScoreUtilities.findAccessionInObo(scoreAccession);
		}
		
		if (cvAccession != null) {
			this.comparator = arg;
			this.value = value;
			this.negate = negate;
			
			this.name = (cvScore != null ? cvScore.getName() : null) + " Filter for PSM";
			this.filteringName = (cvScore != null ? cvScore.getName() : null) + " (PSM)";
			this.shortName = prefix + cvScore.getAccession();
		} else {
			this.comparator = null;
			this.value = null;
			
			this.name = null;
			this.filteringName = null;
			this.shortName = null;
		}
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
			return ((IntermediatePeptideSpectrumMatch) o).getScore(cvAccession);
		} else if (o instanceof Number) {
			return ((Number) o).doubleValue();
		} else {
			return null;
		}
	}
	
	
	@Override
	public boolean supportsClass(Object c) {
		// only support SpectrumIdentification and Number
		if (c instanceof IntermediatePeptideSpectrumMatch) {
			return true;
		} else if (c instanceof Number) {
			return true;
		}
		
		return false;
	}
	
	
	@Override
	public String toString() {
		StringBuilder str = new StringBuilder(getShortName());
		
		if (cvAccession != null) {
			str.append(" (");
			str.append(CvScore.getCvRefByAccession(cvAccession).getName());
			str.append(")");
		}
		
		if (getFilterNegate()) {
			str.append(" not");
		}
		
		str.append(" ");
		str.append(getFilterComparator().toString());
		
		str.append(" ");
		str.append(getFilterValue());
		
		return str.toString();
	}
}
