package be.svlandeg.diffany.concepts;

/**
 * An edge definition holds all information of an edge, except its actual source and target nodes.
 * It is thus a virtual definition of an edge, free from any network context.
 * 
 * It is used by the EdgeOntology to reason about edge types etc. without considering the actual nodes.
 * 
 * @author Sofie Van Landeghem
 */
public class EdgeDefinition
{
	
	protected String type;
	
	protected boolean symmetrical;
	
	protected double weight;	
	protected boolean negated;
	
	/**
	 * Create a new edge definition. It will be initalized to void (i.e. weight = 0).
	 */
	public EdgeDefinition()
	{
		this("*NOTYPE*" , true, 0, false);
	}
	
	/**
	 * Create a new edge definition with a certain type, direction, weight and directionality
	 * @param type the interaction type of this edge
	 * @param symmetrical defines whether the edge is symmetrical or directed
	 * @param weight the weight or confidence of this edge (should be positive)
	 * @param negated defines whether or not the edge is negated (e.g. does NOT bind)
	 * @throws IllegalArgumentException when the specified weight is a negative number
	 */
	public EdgeDefinition(String type, boolean symmetrical, double weight, boolean negated) throws IllegalArgumentException
	{
		setType(type);
		makeSymmetrical(symmetrical);
		setWeight(weight);
		makeNegated(negated);
	}
	
	/**
	 * Obtain a void edge, for the purpose of being able to compare it to existing edges.
	 * @return a void edge (weight == 0, symmetrical == true)
	 */
	public static EdgeDefinition getVoidEdge()
	{
		return new EdgeDefinition("*NOTYPE*" , true, 0, false);
	}
	
	/**
	 * Set the type of the edge
	 * @param type the type of the edge
	 */
	public void setType(String type)
	{
		this.type = type;
	}
	
	/**
	 * Get the type of this edge
	 * @return the edge type
	 */
	public String getType()
	{
		return type;
	}
	
	/**
	 * Make this edge symmetrical or not
	 * @param symmetrical whether or not this edge is symmetrical
	 */
	public void makeSymmetrical(boolean symmetrical)
	{
		this.symmetrical = symmetrical;
	}
	
	/**
	 * Return whether or not this edge is symmetrical 
	 * (if not, it goes specifically from source to target, otherwise there is no real direction)
	 * @return whether or not this edge is symmetrical
	 */
	public boolean isSymmetrical()
	{
		return symmetrical;
	}
	
	/**
	 * Make this edge negated or not
	 * @param negated whether or not this edge is negated
	 */
	public void makeNegated(boolean negated)
	{
		this.negated = negated;
	}
	
	/**
	 * Return whether or not this edge is negated (e.g. does NOT bind) 
	 * @return whether or not this edge is negated
	 */
	public boolean isNegated()
	{
		return negated;
	}
	
	/**
	 * Get the weight of this edge.
	 * When it never has been set, it is 1.0 by default.
	 * The edge weight should always be positive. When it is 0, the edge can be ignored.
	 * @return the edge weight
	 */
	public double getWeight()
	{
		return weight;
	}
	
	/**
	 * Set the weight of this edge (only when the value is appropriate).
	 * @param weight the ne weight of the edge
	 * @throws IllegalArgumentException when the specified weight is a negative number
	 */
	public void setWeight(double weight) throws IllegalArgumentException
	{
		checkWeight(weight);
		this.weight = weight;
	}
	
	/**
	 * Internal method to check whether the weight value is appropriate.
	 * @throws IllegalArgumentException when the specified weight is a negative number
	 */
	private void checkWeight(double weight) throws IllegalArgumentException
	{
		if (weight < 0.0)
		{
			String errormsg = "The edge weight should be positive!";
			throw new IllegalArgumentException(errormsg);
		}
	}

}
