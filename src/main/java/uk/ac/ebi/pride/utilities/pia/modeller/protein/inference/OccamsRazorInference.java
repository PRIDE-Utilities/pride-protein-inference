package uk.ac.ebi.pride.utilities.pia.modeller.protein.inference;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateGroup;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateStructure;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.peptide.PeptideScoring;
import uk.ac.ebi.pride.utilities.pia.modeller.scores.protein.ProteinScoring;




/**
 * This inference filter reports all the PIA {@link IntermediateGroup}s as protein,
 * which together fulfill the Occam's Razor constraints. I.e. the minimal set of
 * groups is reported, which contain all peptides and each protein is contained,
 * at least as subset of peptides.
 * <p>
 * There are 3 constraints for a Group, to be reported:<br/>
 * 1) it has no parent-Groups (but implicit accessions, there cannot be a Group
 * without parents and without accessions)<br/>
 * 2) Group fulfills 1) and has any direct peptides -> report it<br/>
 * 3) Group fulfills 1) but not 2): get the peptide children groups (i.e. the
 * Groups, which have the peptides). If the pepChildGroups are not fully
 * explained by any other Group fulfilling 1), report the Group. (If it is
 * explained by any other, set it as a subGroup of it).
 * 
 * @author julian
 * 
 */
public class OccamsRazorInference extends AbstractProteinInference {
	
	/** the logger for this class */
	private static final Logger logger=  LoggerFactory.getLogger(OccamsRazorInference.class);
	
	/** the human readable name of this filter */
	protected static final String name = "Occam's Razor";
	
	/** the machine readable name of the filter */
	protected static final String shortName = "inference_occams_razor";
	
	/** the progress of the inference */
	private Double progress;
	
	/** variable to calculate the progress faster*/
	private Double progressStep;
	
	
	/** this iterator iterates over the mapping from the tree ID to its groups*/
	private Iterator<Set<IntermediateGroup>> clustersIterator;
	
	/** this list holds the reported proteins */
	private List<InferenceProteinGroup> reportProteinGroups;
	
	
	public OccamsRazorInference(IntermediateStructure intermediateStructure,
			PeptideScoring peptideScoring, ProteinScoring proteinScoring,
			List<AbstractFilter> filters, Integer nrThreads) {
		super(intermediateStructure, peptideScoring, proteinScoring, filters, nrThreads);
		
		this.progress = 0.0;
	}
	
	/*
	@Override
	public List<LabelValueContainer<String>> getFilterTypes() {
		List<LabelValueContainer<String>> filters = new ArrayList<LabelValueContainer<String>>();
		
		// PSM filters
		filters.add(new LabelValueContainer<String>(null, "--- PSM ---"));
		for (Map.Entry<String, String>  scoreIt
				: getAvailableScoreShorts().entrySet()) {
			String[] filterNames = PSMScoreFilter.getShortAndFilteringName(
					scoreIt.getKey(), scoreIt.getValue());
			if (filterNames != null) {
				filters.add(new LabelValueContainer<String>(
						filterNames[0], filterNames[1]));
			}
		}
		filters.add(new LabelValueContainer<String>(NrPSMsPerPSMSetFilter.shortName(),
				NrPSMsPerPSMSetFilter.filteringName()));
		filters.add(new LabelValueContainer<String>(PSMUniqueFilter.shortName(),
				PSMUniqueFilter.filteringName()));
		filters.add(new LabelValueContainer<String>(PSMAccessionsFilter.shortName(),
				PSMAccessionsFilter.filteringName()));
		filters.add(new LabelValueContainer<String>(NrAccessionsPerPSMFilter.shortName(),
				NrAccessionsPerPSMFilter.filteringName()));
		filters.add(new LabelValueContainer<String>(PSMFileListFilter.shortName(),
				PSMFileListFilter.filteringName()));
		filters.add(new LabelValueContainer<String>(PSMModificationsFilter.shortName(),
				PSMModificationsFilter.filteringName()));
		
		//TODO: accessions filter testen!!!!
		
		
		// peptide filters
		filters.add(new LabelValueContainer<String>(null, "--- Peptide ---"));
		filters.add(new LabelValueContainer<String>(NrPSMsPerPeptideFilter.shortName(),
				NrPSMsPerPeptideFilter.filteringName()));
		filters.add(new LabelValueContainer<String>(NrSpectraPerPeptideFilter.shortName(),
				NrSpectraPerPeptideFilter.filteringName()));
		filters.add(new LabelValueContainer<String>(PeptideFileListFilter.shortName(),
				PeptideFileListFilter.filteringName()));
		
		for (Map.Entry<String, String> scoreIt
				: getAvailableScoreShorts().entrySet()) {
			String[] filterNames = PeptideScoreFilter.getShortAndFilteringName(
					scoreIt.getKey(), scoreIt.getValue());
			if (filterNames != null) {
				filters.add(new LabelValueContainer<String>(
						filterNames[0], filterNames[1]));
			}
		}
		
		// protein filters
		filters.add(new LabelValueContainer<String>(null, "--- Protein ---"));
		filters.add(new LabelValueContainer<String>(NrPeptidesPerProteinFilter.shortName(),
				NrPeptidesPerProteinFilter.filteringName()));
		filters.add(new LabelValueContainer<String>(NrPSMsPerProteinFilter.shortName(),
				NrPSMsPerProteinFilter.filteringName()));
		filters.add(new LabelValueContainer<String>(NrSpectraPerProteinFilter.shortName(),
				NrSpectraPerProteinFilter.filteringName()));
		filters.add(new LabelValueContainer<String>(ProteinFileListFilter.shortName(),
				ProteinFileListFilter.filteringName()));
		filters.add(new LabelValueContainer<String>(ProteinScoreFilter.shortName,
				ProteinScoreFilter.filteringName));
		filters.add(new LabelValueContainer<String>(NrUniquePeptidesPerProteinFilter.shortName(),
				NrUniquePeptidesPerProteinFilter.filteringName()));
		
		return filters;
		return null;
	}
	 */
	
	@Override
	public List<InferenceProteinGroup> calculateInference( boolean considerModifications) {
		
		progress = 0.0;
		logger.info(name + " calculateInference started...");
		
		// initialize the cluster iterator
		clustersIterator = intermediateStructure.getClusters().values().iterator();
		
		// initialize the reported list
		reportProteinGroups = new ArrayList<InferenceProteinGroup>();
		
		// the number of threads used for the inference
		List<OccamsRazorWorkerThread> threads =
				new ArrayList<OccamsRazorWorkerThread>(allowedThreads);
		logger.info("using " + allowedThreads + " threads for inference");
		
		// initialize and start  the worker threads
		progressStep = 89.0 / intermediateStructure.getNrClusters();
		for (int i=0; i < allowedThreads; i++) {
			OccamsRazorWorkerThread workerThread =
								new OccamsRazorWorkerThread(i+1,
										this,
										filters,
										considerModifications);
			
			threads.add(workerThread);
			workerThread.start();
		}
		
		progress += 10.0;
		
		// wait for the threads to finish
		for (OccamsRazorWorkerThread workerThread : threads) {
			try {
				workerThread.join();
			} catch (InterruptedException e) {
				// TODO: make better error/exception
				logger.error("thread got interrupted!");
				e.printStackTrace();
			}
		}
		
		progress = 100.0;
		logger.info(name + " calculateInference done");
		return reportProteinGroups;
	}
	
	
	/**
	 * Returns the next cluster in the intermediate structure null, if no more
	 * clusters are available.
	 * 
	 * @return
	 */
	public synchronized Set<IntermediateGroup> getNextCluster() {
		if (clustersIterator != null) {
			if (clustersIterator.hasNext()) {
				return clustersIterator.next();
			}
		}
		return null;
	}
	
	
	/**
	 * Adds the proteinGroups to the list of reported {@link ProteinGroup}s.
	 * 
	 * @param newProteins
	 */
	public synchronized void addToReports(List<InferenceProteinGroup> proteinGroups) {
		reportProteinGroups.addAll(proteinGroups);
		progress += progressStep;
	}
	
	
	@Override
	public String getName() {
		return name;
	}
	
	
	@Override
	public String getShortName() {
		return shortName;
	}
	
	
	@Override
	public Long getProgressValue() {
		return progress.longValue();
	}
}
