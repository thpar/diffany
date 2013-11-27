package be.svlandeg.diffany.semantics;

import java.util.Set;

import be.svlandeg.diffany.concepts.EdgeDefinition;

/**
 * This edge ontology deals with process edges such as ppi and ptm, their interrelationships and their corresponding weights.
 * 
 * @author Sofie Van Landeghem
 */
public class ProcessEdgeOntology extends EdgeOntology
{
	
	/**
	 * Create a new ontology, defining the set of categories. and inserting
	 * After the constructor is called, default edge-category mappings should be inserted using addCategoryMapping!
	 */
	public ProcessEdgeOntology(Set<String> cats)
	{
		removeAllCategoriesAndMappings();
		addCategories(cats);
	}


	@Override
	public EdgeDefinition getDifferentialEdge(EdgeDefinition referenceEdge, EdgeDefinition conditionEdge, double cutoff) throws IllegalArgumentException
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public EdgeDefinition getSharedEdge(EdgeDefinition referenceEdge, EdgeDefinition conditionEdge, double cutoff) throws IllegalArgumentException
	{
		// TODO Auto-generated method stub
		return null;
	}

}
