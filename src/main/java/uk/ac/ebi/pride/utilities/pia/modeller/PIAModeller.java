package uk.ac.ebi.pride.utilities.pia.modeller;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import uk.ac.ebi.pride.utilities.pia.intermediate.DataImportController;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateStructure;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateStructureCreator;
import uk.ac.ebi.pride.utilities.pia.intermediate.prideimpl.PrideImportController;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.psm.PSMModeller;


/**
 * This modeller handles all PIA data, started by the import of identifications,
 * intermediate structure creation and displaying of PSMs, peptides and
 * proteins.
 *  
 * @author julian
 *
 */
public class PIAModeller {
	
	/** logger for this class */
	private static final Logger logger = Logger.getLogger(PIAModeller.class);
	
	
	/** maps from an internal fileID to the corresponding {@link DataImportController} */
	private Map<Integer, DataImportController> inputControllers;
	
	
	/** the allowed number of threads run by PIA */
	int allowedThreads = 4;
	
	
	/** the PSM modeller */
	private PSMModeller psmModeller;
	
	
	/** the intermediate structure creator, not used when loading from intermediate file */
	private IntermediateStructureCreator structCreator;
	
	/** the intermediate structure, either loaded or created */
	private IntermediateStructure intermediateStructure;
	
	
	/**
	 * Creates a modeller for the new creation of a intermediate structure
	 */
	public PIAModeller() {
		inputControllers = new HashMap<Integer, DataImportController>();
		psmModeller = null;
		
		structCreator = new IntermediateStructureCreator(allowedThreads);
		intermediateStructure = null;
	}
	
	
	/**
	 * Creates a modeller which loads the intermediate structure from a file
	 */
	public PIAModeller(String pathname) {
		// the struct creator is not needed for a loaded file
		structCreator = null;
		
		// TODO: load the structure
		inputControllers = null;
		psmModeller = null;
		intermediateStructure = null;
	}
	
	
	/**
	 * As some controllers should be closed, do this here.
	 */
	public void close() {
		// close the inputControllers
		for (DataImportController controller : inputControllers.values()) {
			controller.close();
		}
	}
	
	
	/**
	 * Adds a file to the input files and adds the filtered PSMs to the
	 * structure creator. Filtering is ok, if the used inference methods
	 * are not interfered by it.
	 * <p>
	 * This must be called before the intermediate structure is created.
	 * 
	 * @param pathname
	 * @return
	 */
	public DataImportController addFile(String pathname, List<AbstractFilter> filters) {
		if (structCreator == null) {
			logger.error("the intermediate structure is already created, no more files can be added");
			return null;
		}
		
		File inputFile = new File(pathname);
		
		logger.debug("adding " + inputFile.getAbsolutePath() + " to files");
		
		DataImportController importController =
				new PrideImportController(inputFile, filters);
		// TODO: add the import from other file types and controllers
		
		logger.info("start importing data from the controller ----");
        importController.addAllSpectrumIdentificationsToStructCreator(structCreator);
        
        inputControllers.put(inputControllers.size()+1, importController);
		return importController;
	}
	
	
	/**
	 * Adds a file to the input files and adds all the PSMs to the structure
	 * creator.
	 * <p>
	 * This must be called before the intermediate structure is created.
	 * 
	 * @param pathname
	 * @return
	 */
	public DataImportController addFile(String pathname) {
		return addFile(pathname, null);
	}
	
	
	/**
	 * This method builds the intermediate structure with the data of the input
	 * files. The function can be called only once for a {@link #PIAModeller()},
	 * after all files are added.
	 * 
	 * @return
	 */
	public IntermediateStructure buildIntermediateStructure() {
		logger.debug("starting buildIntermediateStructure");
		
		if (intermediateStructure != null) {
			logger.warn("There is already an intermediate structure created!");
		}
		
		if (structCreator == null) {
			logger.error("The intermediate structure cannot be created, the creator is null!");
			return null;
		}
		
		intermediateStructure =
				structCreator.buildIntermediateStructure();
        
		structCreator = null;
		
		// initialize the PSM modeller
		logger.debug("initializing PSM modeller");
		initializePSMModeller();
		
		logger.debug("buildIntermediateStructure done");
		return intermediateStructure;
	}
	
	
	/**
	 * This method initializes the PSM modeller with the PSMs. The method must
	 * be called after the intermediate structure is built or loaded from file.
	 */
	private void initializePSMModeller() {
		psmModeller = new PSMModeller(inputControllers.size());
		
		// get a mapping from the controllerIDs to the internal fileIDs
		Map<Comparable, Integer> controllerIDtoFileID =
				new HashMap<Comparable, Integer>(inputControllers.size());
		for (Map.Entry<Integer, DataImportController> controllerIt
				: inputControllers.entrySet()) {
			controllerIDtoFileID.put(controllerIt.getValue().getID(),
					controllerIt.getKey());
		}
		
		// distribute the PSMs to the files in the modeller
		for (IntermediatePeptideSpectrumMatch iPSM
				: intermediateStructure.getAllIntermediatePSMs()) {
			Integer fileID = controllerIDtoFileID.get(iPSM.getControllerID());
			psmModeller.addPSMforFile(fileID, iPSM);
		}
		
		
		for (Map.Entry<Integer, DataImportController> controllerIt
				: inputControllers.entrySet()) {
			logger.debug("#PSMs of " + controllerIt.getValue().getInputFileName() +
					": " + psmModeller.getNrPSMs(controllerIt.getKey()) +
					", main score: " + psmModeller.getFdrScoreAccession(controllerIt.getKey()));
		}
	}
}
