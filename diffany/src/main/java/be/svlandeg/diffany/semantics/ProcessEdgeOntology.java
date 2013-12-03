package be.svlandeg.diffany.semantics;

import java.util.Set;

import be.svlandeg.diffany.concepts.EdgeDefinition;

/**
 * This edge ontology deals with process edges such as ppi and ptm, their
 * interrelationships and their corresponding weights.
 * 
 * @author Sofie Van Landeghem
 */
public class ProcessEdgeOntology extends EdgeOntology
{

	private String negPrefix;
	private String posPrefix;

	/**
	 * Create a new ontology, defining the set of categories. and inserting
	 * After the constructor is called, default edge-category mappings should be
	 * inserted using addCategoryMapping!
	 */
	public ProcessEdgeOntology(String posPrefix, String negPrefix, Set<String> cats)
	{
		super();
		this.negPrefix = negPrefix;
		this.posPrefix = posPrefix;
		addCategories(cats);
	}

	@Override
	public EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, EdgeDefinition conEdge, double cutoff) throws IllegalArgumentException
	{
		EdgeDefinition diff_edge = new EdgeDefinition();
		
		boolean refNeg = refEdge.isNegated();
		boolean conNeg = conEdge.isNegated();
		
		if (refNeg)
			refEdge = EdgeDefinition.getVoidEdge();
		
		if (conNeg)
			conEdge = EdgeDefinition.getVoidEdge();

		String refCat = getCategory(refEdge.getType());
		String conCat = getCategory(conEdge.getType());

		boolean equalCats = refCat.equals(conCat);
		Boolean up = null;
		double diffWeight = 0;

		if (equalCats) 
		{
			diffWeight = conEdge.getWeight() - refEdge.getWeight();
			up = true; // assume an increase
			if (diffWeight < 0) // it was a decrease!
			{
				diffWeight *= -1;
				up = false;
			}
		}
		
		String baseType = conCat;

		if (refCat.equals(VOID_TYPE) && ! conCat.equals(VOID_TYPE))
		{
			up = true;
			diffWeight = conEdge.getWeight();
		}

		if (conCat.equals(VOID_TYPE) && ! refCat.equals(VOID_TYPE))
		{
			up = false;
			diffWeight = refEdge.getWeight();
			baseType = refCat;
		}

		
		
		String type = "";
		if (up == null)
		{
			type = refCat + "_to_" + conCat;
			diffWeight = conEdge.getWeight() + refEdge.getWeight();
			if (isChildOf(refCat, conCat))	// the conditional edge is less specific, i.e. not actually differential!
			{
				return EdgeDefinition.getVoidEdge();
			}
			if (isChildOf(conCat, refCat))	// the conditional edge is more specific, i.e. not actually differential!
			{
				return EdgeDefinition.getVoidEdge();
			}
		}
		else 
		{
			if (up)
			{
				type = posPrefix;
			} 
			else
			{
				type = negPrefix;
			}
			type += baseType;
		}
		diff_edge.setType(type);
		
		if (diffWeight <= cutoff)
		{
			return EdgeDefinition.getVoidEdge();
		}

		boolean refSymm = refEdge.isSymmetrical();
		boolean conSymm = conEdge.isSymmetrical();
		boolean diffSymm = refSymm && conSymm;

		diff_edge.setWeight(diffWeight);
		diff_edge.makeSymmetrical(diffSymm);
		diff_edge.makeNegated(false);	// a differential edge is never negated 
		return diff_edge;
	}

}
