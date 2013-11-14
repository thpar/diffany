package be.svlandeg.diffany.semantics;

import java.util.HashSet;
import java.util.Set;

/**
 * Create an Edge Ontology that provides ontological mapping by a simple application of tripples that allow translation of the
 * reference + condition-specific edge category to the differential edge category.
 * Categories will always be matched case independently.
 * 
 * @author Sofie Van Landeghem
 */
public class TrippleEdgeOntology extends EdgeOntology
{

	protected Set<EdgeTripple> tripples;

	/**
	 * Create a tripple edge ontology with an empty set of tripples.
	 */
	public TrippleEdgeOntology()
	{
		super();
		tripples = new HashSet<EdgeTripple>();
	}
	
	/**
	 * Add a new translation tripple to this ontology. Categories are stored case independently.
	 * 
	 * @param referenceCat the category of the reference edge
	 * @param conditionCat the category of the condition-specific edge
	 * @param differentialCat the category of the differential edge
	 * @throws IllegalArgumentException if either of the categories is null
	 */
	public void addTripple(String referenceCat, String conditionCat, String differentialCat) throws IllegalArgumentException
	{
		if (referenceCat == null || conditionCat == null || differentialCat == null)
		{
			String errormsg = "None of the categories in the tripple should be null!";
			throw new IllegalArgumentException(errormsg);
		}
		tripples.add(new EdgeTripple(referenceCat.toLowerCase(), conditionCat.toLowerCase(), differentialCat.toLowerCase()));
	}
	
	/**
	 * Remove all translation tripples.
	 * This does NOT remove all type-to-category mappings (use method removeAllCategories for that).
	 */
	public void removeAllTripples()
	{
		tripples = new HashSet<EdgeTripple>();
	}

	@Override
	public String getDifferentialCategory(String referenceCat, String conditionCat) throws IllegalArgumentException
	{
		if (referenceCat == null || conditionCat == null)
		{
			String errormsg = "None of the categories should be null!";
			throw new IllegalArgumentException(errormsg);
		}
		referenceCat = referenceCat.toLowerCase();
		conditionCat = conditionCat.toLowerCase();
		
		if (referenceCat.equals(conditionCat))
		{
			return null;
		}
		for (EdgeTripple t : tripples)
		{
			if (t.getReferenceCat().equals(referenceCat) && t.getConditionCat().equals(conditionCat))
			{
				return t.getDifferentialCat();
			}
		}
		return referenceCat + "_to_" + conditionCat;
	}

}
