package uk.ac.ebi.pride.utilities.pia.intermediate.prideimpl;

import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.pride.utilities.pia.psisupport.CVUtilities;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;


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
}
