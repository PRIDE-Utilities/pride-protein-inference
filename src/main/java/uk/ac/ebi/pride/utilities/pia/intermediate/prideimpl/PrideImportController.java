package uk.ac.ebi.pride.utilities.pia.intermediate.prideimpl;

import java.io.File;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.pride.utilities.pia.intermediate.DataImportController;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptide;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateProtein;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateStructureCreator;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.FilterUtilities;
import uk.ac.ebi.pride.utilities.data.controller.DataAccessController;
import uk.ac.ebi.pride.utilities.data.controller.impl.ControllerImpl.MzIdentMLControllerImpl;
import uk.ac.ebi.pride.utilities.data.controller.impl.ControllerImpl.PrideXmlControllerImpl;


public class PrideImportController implements DataImportController {
	
	/** the logger for this class */
	private static final Logger logger = LoggerFactory.getLogger(PrideImportController.class);
	
	/** the data access controller */
	private DataAccessController controller;
	
	/** the input file name of this controller */
	private String inputFileName;
	
	/** whether this importController opened the controller and thus has to close it */
	private boolean openedController;
	
	
	/**
	 * Creates an import controller for the given file type.
	 * 
	 * @param inputFile
	 * @param fileType
	 */
	public PrideImportController(File inputFile, InputFileType fileType) {
		initialize(inputFile, fileType);
	}
	
	
	/**
	 * Creates an import controller guessing the file type
	 * 
	 * @param inputFile
	 * @param fileType
	 */
	public PrideImportController(File inputFile) {
		if (MzIdentMLControllerImpl.isValidFormat(inputFile)) {
			initialize(inputFile, InputFileType.MZIDENTML);
		} else if (PrideXmlControllerImpl.isValidFormat(inputFile)) {
			initialize(inputFile, InputFileType.PRIDE_XML);
		} else {
			logger.error("This is not a valid file type!");
        }
	}
	
	
	/**
	 * Creates an import controller for the given {@link DataAccessController}
	 * 
	 * @param inputFile
	 * @param fileType
	 */
	public PrideImportController(DataAccessController controller) {
		this.controller = controller;
		inputFileName = controller.getName();
	}
	
	
	/**
	 * Initializes the class
	 * 
	 * @param inputFile
	 * @param fileType
	 * @param filters
	 */
	private void initialize(File inputFile, InputFileType fileType) {
		inputFileName = null;
		
		switch (fileType) {
			case PRIDE_XML:
				this.controller = new PrideXmlControllerImpl(inputFile);
				break;
				
			case MZIDENTML:
			default:
				this.controller = new MzIdentMLControllerImpl(inputFile, false, true);
				break;
		}
		
		if (this.controller != null) {
			openedController = true;
			inputFileName = inputFile.getAbsolutePath();
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
	
	
	/**
	 * Returns the actual controller.
	 * @return
	 */
	public DataAccessController getController() {
		return controller;
	}
	
	
	@Override
	public void close() {
		if (openedController) {
			controller.close();
		}
	}
	
	
	@Override
	public void addSpectrumIdentificationsToStructCreator(IntermediateStructureCreator structCreator, List<AbstractFilter> filters) {
		int nrProteins = controller.getNumberOfProteins();
		
		logger.info("start importing data from the controller, " + nrProteins +
				" proteins to go");
		
		int processedProtIDs = 0;
		for (Comparable proteinId : controller.getProteinIds()) {
			addProteinsSpectrumIdentificationsToStructCreator(proteinId, structCreator, filters);

			processedProtIDs++;
			if (((processedProtIDs % 1000) == 0) && (processedProtIDs > 1)) {
				logger.info("processed proteins " + processedProtIDs + " / " + nrProteins);
			}
		}
	}
	
	
	/**
	 * Adds the spectrum identifications of a single protein to the structure
	 * creator.
	 * 
	 * @param proteinId
	 * @param structCreator
	 * @return the protein accession of the inserted protein
	 */
	public String addProteinsSpectrumIdentificationsToStructCreator(Comparable proteinId, IntermediateStructureCreator structCreator, List<AbstractFilter> filters) {
		// create the protein, add it later (when there is a filtered PSM)
		IntermediateProtein protein = new PrideIntermediateProtein(controller, proteinId);
		
		for (Comparable peptideId : controller.getPeptideIds(proteinId)) {
			// add the peptides
			IntermediatePeptideSpectrumMatch psm =
					new PrideIntermediatePeptideSpectrumMatch(controller, proteinId, peptideId);
			
			if ((filters == null) || (filters.size() == 0) || (FilterUtilities.satisfiesFilterList(psm, filters))) {
				// add the protein (only, if any PSM passes filters)
				String proteinAccession = protein.getAccession();
				if (!structCreator.proteinsContains(proteinAccession)) {
					structCreator.addProtein(protein);
				} else {
					structCreator.addProteinInformation(proteinAccession, protein);
				}
				
				String pepSequence = psm.getSequence();
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
				structCreator.addPeptideToProteinConnection(pepID, proteinAccession);
			}
		}
		
		return protein.getAccession();
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
