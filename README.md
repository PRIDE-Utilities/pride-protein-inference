pride-protein-inference
===============

# About pride-protein-inference

The primary purpose of pride-protein-inference library is to provide different algorithms for protein inference in PRIDE Proteomics Experiments. You may also find it useful for your own computational proteomics projects.

# License

pride-protein-inference is a PRIDE API licensed under [Apache License 2.0](http://www.apache.org/licenses/LICENSE-2.0.txt).

# How to cite it:

 
# Main Features
* Data structure to handle protein inference such as Protein Groups, Peptides, Scores, Protein Inference Algorithms, etc 
* Functions to compute Protein Inference Algorithms. 
* Set of CVterms for Protein Inference Algorithms such as Protein Scores, PSMs scores, etc
* Filtering PSM, Peptides and Proteins. 

# The library provides four key modules:

**Note**: the library is still evolving, we are committed to expand this library and add more useful classes.

# Getting pride-protein-inference

The zip file in the releases section contains the PRIDE Utilities jar file and all other required libraries.

Maven Dependency

PRIDE Utilities library can be used in Maven projects, you can include the following snippets in your Maven pom file.
 
 ```maven
 <dependency>
   <groupId>uk.ac.ebi.pride.utilities</groupId>
   <artifactId>pride-protein-inference</artifactId>
   <version>0.0.1-SNAPSHOT</version>
 </dependency> 
 ```
 ```maven
 <!-- EBI repo -->
 <repository>
     <id>nexus-ebi-repo</id>
     <url>http://www.ebi.ac.uk/intact/maven/nexus/content/repositories/ebi-repo</url>
 </repository>
 
 <!-- EBI SNAPSHOT repo -->
 <snapshotRepository>
    <id>nexus-ebi-repo-snapshots</id>
    <url>http://www.ebi.ac.uk/intact/maven/nexus/content/repositories/ebi-repo-snapshots</url>
 </snapshotRepository>
```
Note: you need to change the version number to the latest version.

For developers, the latest source code is available from our SVN repository.

# Getting Help

If you have questions or need additional help, please contact the PRIDE Helpdesk at the EBI: pride-support at ebi.ac.uk (replace at with @).

Please send us your feedback, including error reports, improvement suggestions, new feature requests and any other things you might want to suggest to the PRIDE team.

# This library has been used in:


How to use pride-protein-inference
===============

# Using pride-protein-inference 

Here we will show you a simple example to compute protein inference in a mzIdentML file:

### Compute Protein Inference:


```java 
// ---------------------------------------------------------------------

int allowedThreads = 4;
boolean considerModifications = false;
boolean filterPSMsOnImport = false;
String fdrScoreAccession = CvScore.PSI_MASCOT_SCORE.getAccession();
String peptideScoreAccession = CvScore.PSI_MASCOT_SCORE.getAccession();
boolean oboLookup = false;

Logger logger = Logger.getLogger("log.txt");

inputFile = new File("sample-file-without-protein-group.mzid");


// some testing variables
logger.info("using " + inputFile.getAbsolutePath());
		
		
// set up some filters
List<AbstractFilter> filters = new ArrayList<AbstractFilter>();
		
AbstractFilter filter = new PSMScoreFilter(FilterComparator.less_equal, 0.01, false, CvScore.PSI_PSM_LEVEL_FDRSCORE.getAccession(), oboLookup);
filters.add(filter);
		
filter = new PSMDecoyFilter(FilterComparator.equal, false, false);
filters.add(filter);
		
		
if (filterPSMsOnImport) {
     importController = new PrideImportController(inputFile, filters);
} else {
     importController = new PrideImportController(inputFile);
}
        
IntermediateStructureCreator structCreator = new IntermediateStructureCreator(allowedThreads);
		
logger.info("start importing data from the controller");
importController.addAllSpectrumIdentificationsToStructCreator(structCreator);
        
        
 logger.info("creating intermediate structure with\n\t"
		+ structCreator.getNrSpectrumIdentifications() + " spectrum identifications\n\t"
		+ structCreator.getNrPeptides() + " peptides\n\t"
		+ structCreator.getNrProteins() + " protein accessions");
		
IntermediateStructure intermediateStructure = structCreator.buildIntermediateStructure();
		
// sort the PSMs by score
		
logger.info("sorting PSMs by score");
List<IntermediatePeptideSpectrumMatch> psms = new ArrayList<IntermediatePeptideSpectrumMatch>(intermediateStructure.getAllIntermediatePSMs());

logger.info("   obtained PSMs for sorting");
		
Collections.sort(psms, new IntermediatePSMComparator(fdrScoreAccession, oboLookup));
logger.info("sorting done");
		

// then calculate the FDR and FDRScore		
FDRUtilities.calculateFDR(psms, fdrScoreAccession);
logger.info("   fdr calculation done");
		
FDRUtilities.calculateFDRScore(psms,  fdrScoreAccession, false);
logger.info("   FDRScore calculation done");
		
		
// perform the protein inference
		
PeptideScoring pepScoring = new PeptideScoringUseBestPSM(peptideScoreAccession, oboLookup);
ProteinScoring protScoring = new ProteinScoringAdditive(false, pepScoring);
		
AbstractProteinInference proteinInference = new OccamsRazorInference(intermediateStructure, pepScoring, protScoring, filters, allowedThreads);
		
List<InferenceProteinGroup> inferenceGroups = proteinInference.calculateInference(considerModifications);
		
logger.info("inferred groups: " + inferenceGroups.size());
        
        
// print out some information

for (InferenceProteinGroup group : inferenceGroups) {
	  StringBuilder groupText = new StringBuilder();
	  for (IntermediateProtein protein : group.getProteins()) {
		    if (groupText.length() > 0) {
				groupText.append(',');
			}
				groupText.append(protein.getAccession());
			}
			Set<String> subAccessions = new HashSet<String>();
			for (InferenceProteinGroup subGroup : group.getSubGroups()) {
				for (IntermediateProtein subProt : subGroup.getProteins()) {
					subAccessions.add(subProt.getAccession());
				}
			}
			for (String subAccession : subAccessions) {
				groupText.append(',');
				groupText.append('[');
				groupText.append(subAccession);
				groupText.append(']');
			}
			groupText.append(" (");
			groupText.append(group.getScore());
			groupText.append(')');
			logger.info(groupText);
	}
 }
importController.close();

```


Tip: In the next future we will provide more example that you might find them useful. 
