package be.svlandeg.diffany.semantics;

/**
 * An EdgeTripple can translate a reference + condition-specific edge into a differential edge,
 * only by looking at the types. Used by TrippleEdgeOntology.
 * 
 * @author Sofie Van Landeghem
 */
public class EdgeTripple
{
	protected String referenceCat;
	protected String conditionCat;
	protected String differentialCat;
	
	/**
	 * Create a new tripple for translation into a differential category
	 * @param referenceCat the category of the reference edge
	 * @param conditionCat the category of the condition-specific edge
	 * @param differentialCat the category of the differential edge
	 */
	public EdgeTripple(String referenceCat, String conditionCat, String differentialCat)
	{
		this.referenceCat = referenceCat;
		this.conditionCat = conditionCat;
		this.differentialCat = differentialCat;
	}
	
	public boolean equals(Object o)
	{
		if (o instanceof EdgeTripple)
		{
			EdgeTripple t2 = (EdgeTripple) o;
			return (t2.referenceCat.equals(referenceCat) && t2.conditionCat.equals(conditionCat) && t2.differentialCat.equals(differentialCat));
		}
		return false;
	}
	
	/**
	 * Get the category of the reference edge
	 * @return the category of the reference edge
	 */
	public String getReferenceCat()
	{
		return referenceCat;
	}
	
	/**
	 * Get the category of the condition-specific edge
	 * @return the category of the condition-specific edge
	 */
	public String getConditionCat()
	{
		return conditionCat;
	}
	
	/**
	 * Get the category of the differential edge
	 * @return the category of the differential edge
	 */
	public String getDifferentialCat()
	{
		return differentialCat;
	}
	
}
