package uk.ac.ebi.pride.utilities.pia.modeller;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.pride.utilities.data.controller.DataAccessController;
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
	private static final Logger logger = LoggerFactory.getLogger(PIAModeller.class);
	
	
	/** maps from an internal fileID to the corresponding {@link DataImportController} */
	private Map<Integer, DataImportController> inputControllers;
	
	
	/** the allowed number of threads run by PIA */
	int allowedThreads;
	
	
	/** the PSM modeller */
	private PSMModeller psmModeller;
	
	
	/** the intermediate structure creator, not used when loading from intermediate file */
	private IntermediateStructureCreator structCreator;
	
	/** the intermediate structure, either loaded or created */
	private IntermediateStructure intermediateStructure;
	
	
	
	/**
	 * Creates a modeller for the new creation of a intermediate structure using
	 * as many threads as processors are available.
	 */
	public PIAModeller() {
		this(Runtime.getRuntime().availableProcessors());
	}
	
	
	/**
	 * Creates a modeller for the new creation of a intermediate structure using
	 * only the specified number of threads
	 */
	public PIAModeller(int nrThreads) {
		this.allowedThreads = nrThreads;
		this.structCreator = new IntermediateStructureCreator(allowedThreads);
		
		this.inputControllers = new HashMap<Integer, DataImportController>();
		
		// this will be initialized later
		this.intermediateStructure = null;
		
		// this will be initialized after the intermediate structure is created
		this.psmModeller = null;
	}
	
	
	/**
	 * Creates a modeller which loads the intermediate structure from a file
	 * using as many threads as processors are available.
	 */
	public PIAModeller(String pathname) {
		this(pathname, Runtime.getRuntime().availableProcessors());
	}
	
	
	/**
	 * Creates a modeller which loads the intermediate structure from a file
	 * using only the specified number of threads
	 */
	public PIAModeller(String pathname, int nrThreads) {
		this(nrThreads);
		
		// the struct creator is not needed for a loaded file
		structCreator = null;
		
		// TODO: load the structure and set everything correctly
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
	 * Adds a file to the input files.
	 * <p>
	 * This must be called before the intermediate structure is created.
	 * 
	 * @param pathname
	 * @return the ID of the file
	 */
	public Integer addFile(String pathname) {
		if (structCreator == null) {
			logger.error("the intermediate structure is already created, no more files can be added");
			return null;
		}
		
		File inputFile = new File(pathname);
		
		logger.debug("adding " + inputFile.getAbsolutePath() + " to files");
		
		DataImportController importController = new PrideImportController(inputFile);
		// TODO: add the import from other file types and controllers
		
		Integer fileID = inputControllers.size()+1;
        inputControllers.put(fileID, importController);
		return fileID;
	}
	
	
	/**
	 * Adds a {@link DataAccessController} as source for import. Does the same
	 * as {@link #addFile(String)} for other files.
	 * 
	 * @param controller
	 * @return
	 */
	public Integer addPrideControllerAsInput(DataAccessController controller) {
		if (structCreator == null) {
			logger.error("the intermediate structure is already created, no more files can be added");
			return null;
		}
		
		logger.debug("adding pride controller \"" + controller.getName() +  "\" to files");
		DataImportController importController = new PrideImportController(controller);
		
		Integer fileID = inputControllers.size()+1;
        inputControllers.put(fileID, importController);
		return fileID;
	}
	
	
	/**
	 * Adds a file to the input files and adds the filtered PSMs to the
	 * structure creator. Filtering is ok, if the used inference methods
	 * are not interfered by it.
	 * <p>
	 * This must be called before the intermediate structure is created.
	 * 
	 * @param pathname
	 * @return the ID of the file
	 */
	public Integer addFileAndImportSpectra(String pathname, List<AbstractFilter> filters) {
		Integer fileID = addFile(pathname);
		
		importDataFromFile(fileID, filters);
        
		return fileID;
	}
	
	
	/**
	 * Imports all spectrum identifications from the file to the structure creator
	 * 
	 * @param fileID
	 */
	public void importAllDataFromFile(Integer fileID) {
		importDataFromFile(fileID, null);
	}
	
	
	/**
	 * Imports the filtered PSMs and adds them to the structure creator.
	 * 
	 * @param fileID
	 */
	public void importDataFromFile(Integer fileID, List<AbstractFilter> filters) {
		if (structCreator == null) {
			logger.error("the intermediate structure is already created, no more files can be added");
			return;
		}
		
		inputControllers.get(fileID).addSpectrumIdentificationsToStructCreator(structCreator, filters);
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
		
		logger.debug("buildIntermediateStructure done, #clusters " + intermediateStructure.getNrClusters());
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
	
	
	/**
	 * Returns the {@link DataImportController} of the given file.
	 * 
	 * @param fileID
	 * @return
	 */
	public DataImportController getImportController(Integer fileID) {
		return inputControllers.get(fileID);
	}
	
	
	/**
	 * Returns the {@link IntermediateStructureCreator} used by the modeller.
	 * 
	 * @return
	 */
	public IntermediateStructureCreator getIntermediateStructureCreator() {
		return structCreator;
	}
	
	
	/**
	 * Returns the PSMModeller
	 * <p>
	 * Is initialized after the intermediate structure is created or loaded
	 * @return
	 */
	public PSMModeller getPSMModeller() {
		return psmModeller;
	}
}
