package uk.ac.ebi.pride.utilities.pia.intermediate.prideimpl;

import java.io.File;
import java.util.ArrayList;
import java.util.List;


import org.apache.log4j.Logger;

import uk.ac.ebi.pride.utilities.data.controller.DataAccessController;
import uk.ac.ebi.pride.utilities.data.controller.impl.ControllerImpl.MzIdentMLControllerImpl;
import uk.ac.ebi.pride.utilities.data.controller.impl.ControllerImpl.PrideXmlControllerImpl;
import uk.ac.ebi.pride.utilities.data.utils.MzIdentMLUtils;
import uk.ac.ebi.pride.utilities.pia.intermediate.DataImportController;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptide;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateProtein;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateStructureCreator;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterUtilities;


public class PrideImportController implements DataImportController {
	
	/** the logger for this class */
	private static final Logger logger= Logger.getLogger(PrideImportController.class);
	
	/** the data access controller */
	private DataAccessController controller;
	
	/** the applied filters during the import */
	private List<AbstractFilter> filters;
	
	/** the input file name of this controller */
	private String inputFileName;
	
	
	/**
	 * Creates an import controller for the given file type without any filters
	 * used for import.
	 * 
	 * @param inputFile
	 * @param fileType
	 */
	public PrideImportController(File inputFile, InputFileType fileType) {
		this(inputFile, fileType, null);
	}
	
	
	/**
	 * Creates an import controller guessing the file type and not filtering
	 * on import
	 * 
	 * @param inputFile
	 */
	public PrideImportController(File inputFile) {
		this(inputFile, (List<AbstractFilter>)null);
	}
	
	
	/**
	 * Creates an import controller guessing the file type using the given
	 * filters on PSM import.
	 * 
	 * @param inputFile
	 */
	public PrideImportController(File inputFile, List<AbstractFilter> filters) {
		Class fileType = MzIdentMLUtils.getFileType(inputFile);
		
		if (fileType == MzIdentMLControllerImpl.class) {
			initialize(inputFile, InputFileType.MZIDENTML, null);
		} else if (fileType == PrideXmlControllerImpl.class) {
			initialize(inputFile, InputFileType.PRIDE_XML, null);
		}
	}
	
	
	/**
	 * Creates an import controller for the given file type using the given
	 * filters on PSM import.
	 * 
	 * @param inputFile
	 * @param fileType
	 */
	public PrideImportController(File inputFile, InputFileType fileType, List<AbstractFilter> filters) {
		initialize(inputFile, fileType, filters);
	}
	
	
	/**
	 * Initializes the class
	 * 
	 * @param inputFile
	 * @param fileType
	 * @param filters
	 */
	private void initialize(File inputFile, InputFileType fileType, List<AbstractFilter> filters) {
		inputFileName = null;
		
		switch (fileType) {
			case PRIDE_XML:
				this.controller = new PrideXmlControllerImpl(inputFile);
				break;
				
			case MZIDENTML:
			default:
				this.controller = new MzIdentMLControllerImpl(inputFile, true, true);
				break;
		}
		
		if (this.controller != null) {
			inputFileName = inputFile.getAbsolutePath();
			this.filters = (filters == null) ? new ArrayList<AbstractFilter>() : filters;
			
			controller.getContentCategories();
		}
	}
	
	
	@Override
	public String getID() {
		return controller.getUid();
	}
	
	
	@Override
	public String getInputFileName() {
		return inputFileName;
	}
	
	
	@Override
	public void close() {
		controller.close();
	}
	
	
	@Override
	public void addAllSpectrumIdentificationsToStructCreator(IntermediateStructureCreator structCreator) {
		int nrProteins = controller.getNumberOfProteins();
		logger.info(nrProteins + " proteins to go");
		
		int processedProtIDs = 0;
		for (Comparable proteinId : controller.getProteinIds()) {
			
			// create the protein, add it later (when there is a filtered PSM)
			IntermediateProtein protein = new PrideIntermediateProtein(controller, proteinId);
			
			for (Comparable peptideId : controller.getPeptideIds(proteinId)) {
				// add the peptides
				IntermediatePeptideSpectrumMatch psm =
						new PrideIntermediatePeptideSpectrumMatch(controller, proteinId, peptideId);
				
				if ((filters.size() == 0) || (FilterUtilities.satisfiesFilterList(psm, filters))) {
					// add the protein (only, if any PSM passes filters)
					if (!structCreator.proteinsContains(protein.getID())) {
						structCreator.addProtein(protein);
					}
					
					String pepSequence = psm.getSpectrumIdentification().getSequence();
					Comparable pepID = IntermediatePeptide.computeID(pepSequence);
					
					IntermediatePeptide peptide;
					if (structCreator.peptidesContains(pepID)) {
						peptide = structCreator.getPeptide(pepID);
					} else {
						peptide = new IntermediatePeptide(pepSequence);
						structCreator.addPeptide(peptide);
					}
					
					// add the PSM to the peptide (if it does not already exist)
					peptide.addPeptideSpectrumMatch(psm);
					
					// connect the peptide and protein
					structCreator.addPeptideToProteinConnection(pepID, protein.getID());
				}
			}
			processedProtIDs++;
			if (((processedProtIDs % 1000) == 0) && (processedProtIDs > 1)) {
				logger.info("processed proteins " + processedProtIDs + " / " + nrProteins);
			}
		}
	}
	
	
	/**
	 * Defines the input file type
	 * @author julian
	 *
	 */
	public enum InputFileType {
		MZIDENTML,
		PRIDE_XML,
	}
}
