package uk.ac.ebi.pride.utilities.pia.intermediate;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;


public class IntermediateGroup {
	
	/** ID of the group */
	private Integer ID;
	
	/** ID of the group's tree */
	private Integer treeID;
	
	/** List of the direct peptides in the group. */
	private Set<IntermediatePeptide> peptides;
	
	/** Children groups of this group, i.e. groups where this group points to. */
	private Set<IntermediateGroup> children;
	
	/** Parents of this group, i.e. groups pointing to this group. */
	private Set<IntermediateGroup> parents;
	
	/** List of the directly attached proteins of this group. */
	private Set<IntermediateProtein> proteins;
	
	
	/**
	 * Basic Constructor, sets all the maps to null and score to NaN.
	 * 
	 * @param id
	 */
	public IntermediateGroup(int id) {
		this.ID = id;
		this.treeID = -1;
		this.peptides = null;
		this.children = null;
		this.parents = null;
		this.proteins = null;
	}
	
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
		
		IntermediateGroup group = (IntermediateGroup)obj;
		
		return (ID == group.ID) &&
				(treeID == group.treeID) && 
				((peptides == null) ? (group.peptides == null) : peptides.equals(group.peptides)) &&
				((children == null) ? (group.children == null) : children.equals(group.children)) &&
				((parents == null) ? (group.parents == null) : parents.equals(group.parents)) &&
				((proteins == null) ? (group.proteins == null) : proteins.equals(group.proteins));
	}
	
	
	@Override
	public int hashCode() {
		int result;
		
		result = ID.hashCode();
		result = 31*result + treeID.hashCode();
        result = 31 * result + ((peptides != null) ? peptides.hashCode() : 0);
        
		// bad to take the childrens' hashcode directly, as the group is referenced as parent
		int midResult = 0;
		if (children != null) {
			for (IntermediateGroup child : children) {
				 midResult = midResult*31 + child.getID().hashCode();
			}
		}
		result = 31*result + midResult;
		
		// wbad to take the parents' hashcode directly, as the group is referenced as child
		midResult = 0;
		if (parents != null) {
			for (IntermediateGroup parent : parents) {
				midResult = midResult*31 + parent.getID().hashCode();
			}
		}
		result = 31*result + midResult;
		
        result = 31 * result + ((proteins != null) ? proteins.hashCode() : 0);
        
        return result;
	}
	
	
	/**
	 * Setter for the ID.
	 * 
	 * @return
	 */
	public void setID(Integer id) {
		this.ID = id;
	}
	
	
	/**
	 * Getter for the ID.
	 * 
	 * @return
	 */
	public Integer getID() {
		return ID;
	}
	
	
	/**
	 * Setter for the treeID.
	 * 
	 * @param id
	 */
	public void setTreeID(Integer id) {
		this.treeID = id;
	}
	
	
	/**
	 * Getter for the group's treeID.
	 * 
	 * @return
	 */
	public Integer getTreeID() {
		return treeID;
	}
	
	
	/**
	 * Setter for the peptides.
	 * 
	 * @param peptides
	 */
	public void setPeptides(Collection<IntermediatePeptide> peptides) {
		this.peptides = new HashSet<IntermediatePeptide>(peptides);
	}
	
	
	/**
	 * Adds a single peptide to the group.
	 * 
	 * @param peptides
	 */
	public void addPeptide(IntermediatePeptide peptide) {
		if (peptides == null) {
			peptides = new HashSet<IntermediatePeptide>();
		}
		
		peptides.add(peptide);
	}
	
	
	/**
	 * Getter for the peptides.
	 * 
	 * @return
	 */
	public Set<IntermediatePeptide> getPeptides() {
		return peptides;
	}
	
	
	/**
	 * getter for all peptides, including all children's peptides.
	 * @return
	 */
	public Set<IntermediatePeptide> getAllPeptides() {
		Set<IntermediatePeptide> pepSet = new HashSet<IntermediatePeptide>();
		
		if (peptides != null) {
			for (IntermediatePeptide peptide : peptides) {
				pepSet.add(peptide);
			}
		}
		
		if (children != null) {
			for (IntermediateGroup child : getAllPeptideChildren()) {
				pepSet.addAll(child.peptides);
			}
		}
		
		return pepSet;
	}
	
	
	/**
	 * Adds a child to the set of children.
	 * If the set is not yet initialized, initialize it.
	 * 
	 * @param peptides
	 */
	public void addChild(IntermediateGroup child) {
		if (children == null) {
			children = new HashSet<IntermediateGroup>();
		}
		children.add(child);
	}
	
	
	/**
	 * Getter for the children set
	 * 
	 * @return
	 */
	public Set<IntermediateGroup> getChildren() {
		return children;
	}
	
	
	/**
	 * Getter for all children groups of this group, including children's
	 * children and so on.
	 */
	public Set<IntermediateGroup> getAllChildren(){
		Set<IntermediateGroup> allChildren = new HashSet<IntermediateGroup>();
		
		if (children != null) {
			for (IntermediateGroup child : children) {
				allChildren.add(child);
				
				for (IntermediateGroup childChildren : child.getAllChildren()) {
					allChildren.add(childChildren);
				}
			}
		}
		
		return allChildren;
	}
	
	
	/**
	 * Getter for all children groups of this group that have at least one
	 * peptide, recursive, i.e. get the reporting peptide groups.
	 */
	public Set<IntermediateGroup> getAllPeptideChildren(){
		Set<IntermediateGroup> allChildren = new HashSet<IntermediateGroup>();
		
		if (children != null) {
			for (IntermediateGroup child : children) {
				if ((child.peptides != null) && (child.peptides.size() > 0)) {
					allChildren.add(child);
				}
				
				allChildren.addAll(child.getAllPeptideChildren());
			}
		}
		
		return allChildren;
	}
	
	
	/**
	 * Adds a new group to the set of parents.
	 * If the map is not yet initialized, initialize it.
	 * 
	 * @param parent
	 */
	public void addParent(IntermediateGroup parent) {
		if (parents == null) {
			parents = new HashSet<IntermediateGroup>();
		}
		parents.add(parent);
	}
	
	
	/**
	 * Getter for the parents.
	 * 
	 * @return
	 */
	public Set<IntermediateGroup> getParents() {
		return parents;
	}
	
	
	/**
	 * Adds a new protein to the map of proteins.
	 * <p>
	 * If the map is not yet initialized, initialize it.
	 * 
	 * @param dbSeq
	 */
	public void addProtein(IntermediateProtein protein) {
		if (proteins == null) {
			proteins = new HashSet<IntermediateProtein>();
		}
		
		proteins.add(protein);
	}
	
	
	/**
	 * Removes the given dbSequence from the directly connected dbSequences
	 * 
	 * @param dbSeq
	 */
	public void removeProtein(IntermediateProtein protein) {
		if (proteins != null) {
			proteins.remove(protein);
		}
	}
	
	
	/**
	 * Getter for the proteins.
	 * 
	 * @return
	 */
	public Set<IntermediateProtein> getProteins() {
		return proteins;
	}
	
	
	/**
	 * Getter for the intermediateProtein of this group and all the parents.
	 * 
	 * @return
	 */
	public Set<IntermediateProtein> getAllProteins() {
		Set<IntermediateProtein> allProteins = new HashSet<IntermediateProtein>();
		
		if (parents != null) {
			for (IntermediateGroup parent : parents) {
				allProteins.addAll(parent.getAllProteins());
			}
		}
		
		if (proteins != null) {
			allProteins.addAll(proteins);
		}
		
		return allProteins;
	}
	
	
	/**
	 * Returns true, if the group is directly connected to the dbSequence.
	 * 
	 * @return
	 */
	public boolean directlyConnectedToProtein(IntermediateProtein protein) {
		if (proteins != null) {
			return proteins.contains(protein);
		} else {
			return false;
		}
	}
}
