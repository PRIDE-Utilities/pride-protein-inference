package uk.ac.ebi.pride.utilities.pia.modeller;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.pride.utilities.data.controller.DataAccessController;
import uk.ac.ebi.pride.utilities.pia.intermediate.DataImportController;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediatePeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateStructure;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateStructureCreator;
import uk.ac.ebi.pride.utilities.pia.intermediate.prideimpl.PrideImportController;
import uk.ac.ebi.pride.utilities.pia.modeller.filter.AbstractFilter;
import uk.ac.ebi.pride.utilities.pia.modeller.protein.ProteinModeller;
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
	private int allowedThreads;
	
	/** whether to look for unknown CVs in the online OBO */
	private boolean oboLookup;
	
	
	/** the PSM modeller */
	private PSMModeller psmModeller;
	
	/** the protein modeller */
	private ProteinModeller proteinModeller;
	
	
	/** the intermediate structure creator, not used when loading from intermediate file */
	private IntermediateStructureCreator structCreator;
	
	/** the intermediate structure, either loaded or created */
	private IntermediateStructure intermediateStructure;
	
	
	
	/**
	 * Creates a modeller for the new creation of a intermediate structure using
	 * as many threads as processors are available.
	 */
	public PIAModeller() {
		this(Runtime.getRuntime().availableProcessors(), false);
	}
	
	
	/**
	 * Creates a modeller for the new creation of a intermediate structure using
	 * only the specified number of threads
	 */
	public PIAModeller(int nrThreads, boolean oboLookup) {
		this.allowedThreads = nrThreads;
		this.structCreator = new IntermediateStructureCreator(allowedThreads);
		
		this.inputControllers = new HashMap<Integer, DataImportController>();
		
		// this will be initialized later
		this.intermediateStructure = null;
		
		// these will be initialized after the intermediate structure is created
		this.psmModeller = null;
		this.proteinModeller = null;
		
		this.oboLookup = oboLookup;
	}
	
	
	/**
	 * Creates a modeller which loads the intermediate structure from a file
	 * using as many threads as processors are available.
	 */
	public PIAModeller(String pathname) {
		this(pathname, Runtime.getRuntime().availableProcessors(), false);
	}
	
	
	/**
	 * Creates a modeller which loads the intermediate structure from a file
	 * using only the specified number of threads
	 */
	public PIAModeller(String pathname, int nrThreads, boolean oboLookup) {
		this(nrThreads, oboLookup);
		
		// the struct creator is not needed for a loaded file
		this.structCreator = null;
		
		// TODO: load the structure and set everything correctly
		this.inputControllers = null;
		this.psmModeller = null;
		this.proteinModeller = null;
		this.intermediateStructure = null;
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
		
		DataImportController importController = new PrideImportController(inputFile, structCreator);
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
		DataImportController importController = new PrideImportController(controller, structCreator);
		
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
		return addFileAndImportSpectra(pathname, filters, null);
	}
	
	
	/**
	 * Adds a file to the input files and adds the filtered PSMs to the
	 * structure creator. Filtering is ok, if the used inference methods
	 * are not interfered by it.
	 * <p>
	 * Additionally, for the given score accession the scores will be cached
	 * for fast access, if the controller supports it.
	 * <p>
	 * This must be called before the intermediate structure is created.
	 * 
	 * @param pathname
	 * @return the ID of the file
	 */
	public Integer addFileAndImportSpectra(String pathname, List<AbstractFilter> filters,
			String cacheScoreAccession) {
		Integer fileID = addFile(pathname);
		importDataFromFile(fileID, filters, cacheScoreAccession);
		return fileID;
	}
	
	
	/**
	 * Imports all spectrum identifications from the file to the structure
	 * creator.
	 * 
	 * @param fileID
	 */
	public void importAllDataFromFile(Integer fileID) {
		importDataFromFile(fileID, null);
	}
	
	
	/**
	 * Imports all spectrum identifications from the file to the structure
	 * creator. If score caching is supported, cache the given score name.
	 * 
	 * @param fileID
	 */
	public void importAllDataFromFile(Integer fileID, String cacheScoreAccession) {
		importDataFromFile(fileID, null, cacheScoreAccession);
	}
	
	
	/**
	 * Imports the filtered PSMs and adds them to the structure creator.
	 * 
	 * @param fileID
	 */
	public void importDataFromFile(Integer fileID, List<AbstractFilter> filters) {
		importDataFromFile(fileID, filters, null);
	}
	
	
	/**
	 * Imports the filtered PSMs and adds them to the structure creator. If
	 * score caching is supported, cache the given score name.
	 * 
	 * @param fileID
	 */
	public void importDataFromFile(Integer fileID, List<AbstractFilter> filters, String cacheScoreAccession) {
		if (structCreator == null) {
			logger.error("the intermediate structure is already created, no more files can be added");
			return;
		}
		
		inputControllers.get(fileID).addSpectrumIdentificationsToStructCreator(filters, cacheScoreAccession);
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
		
		// initialize the protein modeller
		proteinModeller = new ProteinModeller(intermediateStructure, allowedThreads);
		
		logger.debug("buildIntermediateStructure done, #clusters " + intermediateStructure.getNrClusters());
		return intermediateStructure;
	}
	
	
	/**
	 * This method initializes the PSM modeller with the PSMs. The method must
	 * be called after the intermediate structure is built or loaded from file.
	 */
	private void initializePSMModeller() {
		psmModeller = new PSMModeller(inputControllers.size(), oboLookup);
		
		// get a mapping from the controllerIDs to the internal fileIDs
		Map<Comparable, Integer> controllerIDtoFileID =
				new HashMap<Comparable, Integer>(inputControllers.size());
		for (Map.Entry<Integer, DataImportController> controllerIt
				: inputControllers.entrySet()) {
			controllerIDtoFileID.put(controllerIt.getValue().getID(),
					controllerIt.getKey());
		}
		
		// distribute the PSMs to the files in the modeller
		if (controllerIDtoFileID.size() > 1) {
			for (IntermediatePeptideSpectrumMatch iPSM : intermediateStructure.getAllIntermediatePSMs()) {
				Integer fileID = controllerIDtoFileID.get(iPSM.getControllerID());
				psmModeller.addPSMforFile(fileID, iPSM);
			}
		} else {
			Integer fileID = null;
			for (IntermediatePeptideSpectrumMatch iPSM : intermediateStructure.getAllIntermediatePSMs()) {
				if (fileID == null) {
					fileID = controllerIDtoFileID.get(iPSM.getControllerID());
				}
				psmModeller.addPSMforFile(fileID, iPSM);
			}
		}
		
		
		for (Map.Entry<Integer, DataImportController> controllerIt
				: inputControllers.entrySet()) {
			logger.debug("#PSMs of " + controllerIt.getValue().getInputFileName() +
					": " + psmModeller.getNrPSMs(controllerIt.getKey()));
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
	 * Returns the IDs of the input controllers
	 * @return
	 */
	public Set<Integer> getImportControllerIds() {
		return inputControllers.keySet();
	}
	
	
	/**
	 * Returns the {@link IntermediateStructure} created or loaded by this
	 * modeller.
	 * 
	 * @return
	 */
	public IntermediateStructure getIntermediateStructure() {
		return intermediateStructure;
	}
	
	
	/**
	 * Returns the associated PSMModeller
	 * <p>
	 * Is initialized after the intermediate structure is created or loaded
	 * @return
	 */
	public PSMModeller getPSMModeller() {
		return psmModeller;
	}
	
	
	/**
	 * Returns the associated ProteinModeller
	 * <p>
	 * Is initialized after the intermediate structure is created or loaded
	 * @return
	 */
	public ProteinModeller getProteinModeller() {
		return proteinModeller;
	}
}
