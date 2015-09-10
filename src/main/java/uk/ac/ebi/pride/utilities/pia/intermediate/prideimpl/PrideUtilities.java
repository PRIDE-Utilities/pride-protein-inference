package uk.ac.ebi.pride.utilities.pia.intermediate.prideimpl;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.pride.utilities.data.core.Protein;
import uk.ac.ebi.pride.utilities.data.core.ProteinGroup;
import uk.ac.ebi.pride.utilities.pia.intermediate.IntermediateProtein;
import uk.ac.ebi.pride.utilities.pia.modeller.protein.inference.InferenceProteinGroup;
import uk.ac.ebi.pride.utilities.pia.psisupport.CVUtilities;


/**
 * This class contains some helper functions for the pride import
 * @author julian
 *
 */
public class PrideUtilities {
	
	/**
	 * do not instantiate!
	 */
	private PrideUtilities() {
	}
	
	
	/**
	 * Converts a PRIDE-CvParam into a jMzIdentML CvParam
	 * 
	 * @param prideParams
	 * @return
	 */
	public static List<CvParam> convertCvParams(List<uk.ac.ebi.pride.utilities.data.core.CvParam> prideParams) {
		List<CvParam> paramList = new ArrayList<CvParam>(prideParams.size());
		
		for (uk.ac.ebi.pride.utilities.data.core.CvParam prideParam : prideParams) {
			CvParam cvParam = new CvParam();
			
			cvParam.setAccession(prideParam.getAccession());
			cvParam.setCv(CVUtilities.getByRepresentingName(prideParam.getUnitCVLookupID()));
			cvParam.setName(prideParam.getName());
			cvParam.setUnitAccession(prideParam.getUnitAcc());
			cvParam.setUnitCv(CVUtilities.getByRepresentingName(prideParam.getUnitCVLookupID()));
			cvParam.setUnitName(prideParam.getUnitName());
			cvParam.setValue(prideParam.getValue());
			
			paramList.add(cvParam);
		}
		
		return paramList;
	}
	
	
	/**
	 * Converts a PRIDE-UserParam into a jMzIdentML-UserParam
	 * 
	 * @param prideParams
	 * @return
	 */
	public static List<UserParam> convertUserParams(List<uk.ac.ebi.pride.utilities.data.core.UserParam> prideParams) {
		List<UserParam> paramList = new ArrayList<UserParam>(prideParams.size());
		
		for (uk.ac.ebi.pride.utilities.data.core.UserParam prideParam : prideParams) {
			UserParam userParam = new UserParam();
			
			userParam.setName(prideParam.getName());
			userParam.setType(prideParam.getType());
			userParam.setUnitAccession(prideParam.getUnitAcc());
			userParam.setUnitCv(CVUtilities.getByRepresentingName(prideParam.getUnitCVLookupID()));
			userParam.setUnitName(prideParam.getUnitName());
			userParam.setValue(prideParam.getValue());
			
			paramList.add(userParam);
		}
		
		return paramList;
	}
	
	
	/**
	 * Converts a {@link InferenceProteinGroup} from any protein inference into
	 * a PRIDE {@link ProteinGroup}.
	 * 
	 * @param piaGroup the inferred protein group
	 * @param includingSubGroups whether to add all sub-group entries to the {@link Protein}s of the PRIDE {@link ProteinGroup}
	 * @param completePrideProteins whether to add complete PRIDE proteins (exactly as when no filters were used) or create a new {@link Protein} with the inference information
	 * @return
	 */
	public static ProteinGroup convertInferenceProteinGroup(InferenceProteinGroup piaGroup, boolean includingSubGroups,
			boolean completePrideProteins) {
		Set<Comparable> proteinIDs = new HashSet<Comparable>();
		List<Protein> proteinList = new ArrayList<Protein>();
		
		for (IntermediateProtein intermediateProtein : piaGroup.getProteins()) {
			Protein protein = convertIntermediateProtein((PrideIntermediateProtein)intermediateProtein, completePrideProteins);
			if (protein != null) {
				proteinList.add(protein);
				proteinIDs.add(((PrideIntermediateProtein)intermediateProtein).getPrideProteinID());
			}
		}
		
		if (includingSubGroups) {
			for (InferenceProteinGroup subGroub : piaGroup.getSubGroups()) {
				for (IntermediateProtein intermediateProtein : subGroub.getProteins()) {
					Protein protein = convertIntermediateProtein((PrideIntermediateProtein)intermediateProtein, completePrideProteins);
					if (protein != null) {
						Comparable proteinID = ((PrideIntermediateProtein)intermediateProtein).getPrideProteinID();
						if (!proteinIDs.contains(proteinID)) {
							proteinList.add(protein);
							proteinIDs.add(proteinID);
						}
					}
				}
			}
		}
		
		return new ProteinGroup(piaGroup.getID(), piaGroup.getID(), proteinList);

	}
	
	
	/**
	 * Converts a {@link PrideIntermediateProtein} into a PRIDE {@link Protein}
	 * 
	 * @param intermediateProtein the interemdiate protein
	 * @param completePrideProtein whether to return the complete PRIDE protein or create a new {@link Protein} with the inference information  
	 * @return
	 */
	public static Protein convertIntermediateProtein(PrideIntermediateProtein intermediateProtein, boolean completePrideProtein) {
		Protein protein = null;
		
		if (completePrideProtein) {
			// return the complete PRIDE protein (used when no filters were used in inference)
			protein = intermediateProtein.getPrideProtein();
		}
		//else {
//			// create the inferred protein
//			// TODO: implement!
//            /*
//            Protein protein = new Protein(intermediateProtein.getAccession(),
//                    intermediateProtein.getAccession(),
//                    dbSequence, passThreshold, peptides, score, threshold, sequenceCoverage, gel);
//            */
//		}
		
		return protein;
	}
}
