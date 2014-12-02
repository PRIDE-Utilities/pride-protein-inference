package uk.ac.ebi.pride.utilities.pia.psisupport;

import java.util.Arrays;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.jmzidml.model.mzidml.Cv;

/**
 * Class containing some helpers for the CV managing in the PSI formats
 * 
 * @author julian
 *
 */
public class CVUtilities {
	
	/** the logger for this class */
	private static final Logger logger = LoggerFactory.getLogger(CVUtilities.class);
	
	
	public static enum CVs {
		PSI {
			public Cv getCV() {
				return psiCV;
			}
			
			@Override
			public List<String> getRepresentations() {
				return Arrays.asList(new String[]{
						"PSI",
						"psi",
						"PSI-MS",
						"psi-ms",
				});
			}
		},
		
		UNIMOD {
			public Cv getCV() {
				return unimodCV;
			}
			
			@Override
			public List<String> getRepresentations() {
				return Arrays.asList(new String[]{
						"unimod",
						"UNIMOD",
				});
			}
		},
		
		UO {
			public Cv getCV() {
				return uoCV;
			}
			
			@Override
			public List<String> getRepresentations() {
				return Arrays.asList(new String[]{
						"uo",
						"UO",
						"UNIT-ONTOLOGY",
						"unit-ontology",
				});
			}
		},
		;
		
		/**
		 * represents the PSI-CV
		 */
		public final static Cv psiCV;
		
		
		/**
		 * represents the UNIMOD
		 */
		public final static Cv unimodCV;
		
		
		/**
		 * represents the UnitOntology CV
		 */
		public final static Cv uoCV;
		
		
		static {
			// initialize the static variables
			psiCV = new Cv();
			psiCV.setFullName("Proteomics Standards Initiative Mass Spectrometry Vocabularies");
			psiCV.setId("PSI-MS");
			psiCV.setUri("http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");
			//psiCV.setVersion(value); // changes so often...
			
			unimodCV = new Cv();
			unimodCV.setFullName("UNIMOD");
			unimodCV.setId("UNIMOD");
			unimodCV.setUri("http://www.unimod.org/obo/unimod.obo");
			
			uoCV = new Cv();
			uoCV.setFullName("The Ontology of Units of Measurement");
			uoCV.setId("UO");
			uoCV.setUri("http://unit-ontology.googlecode.com/svn/trunk/unit.obo");
		}
		
		
		/**
		 * return the CV representative
		 * @return
		 */
		abstract public Cv getCV();
		
		
		/**
		 * Returns a list of string representations which may be used by this CV
		 * @return
		 */
		abstract public List<String> getRepresentations();
	}
	
	
	/**
	 * Retrieves the CV which is represented by the given name. If no such CV is
	 * found, returns null.
	 * 
	 * @return the represented CV, if not found null
	 */
	public static Cv getByRepresentingName(String cvName) {
		if (cvName == null) {
			return null;
		}
		
		for (CVs cv : CVs.values()) {
			if (cv.getRepresentations().contains(cvName)) {
				return cv.getCV();
			}
		}
		
		logger.warn("no CV found for " + cvName);
		return null;
	}
	
	
	/**
	 * do not instantiate
	 */
	private CVUtilities() {
		
	}
}
