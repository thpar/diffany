package be.svlandeg.diffany.core.networks;

import be.svlandeg.diffany.core.io.EdgeIO;


/**
 * An edge definition holds all information of an edge, except its actual source and target nodes.
 * It is thus a virtual definition of an edge, free from any network context.
 * In contrast to an edge, an edge definition can still be changed after construction.
 * 
 * It is used by the EdgeOntology to reason about edge types etc. without considering the actual nodes.
 * 
 * @author Sofie Van Landeghem
 */
public class EdgeDefinition
{
	
	public static String DEFAULT_TYPE = "unspecified_connection";
	public static double DEFAULT_WEIGHT = 1.0;
	public static boolean DEFAULT_SYMM = true;
	public static boolean DEFAULT_NEG = false;
	
	protected String type;
	
	protected boolean symmetrical;
	
	protected double weight;	
	protected boolean negated;
	
	/**
	 * Create a new edge definition with a certain type, direction, weight and directionality
	 * 
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
	 * Cloning constructor
	 * @param old the EdgeDefinition to be cloned
	 */
	public EdgeDefinition(EdgeDefinition old)
	{
		this(new String(old.type) , new Boolean(old.symmetrical), old.weight, new Boolean(old.negated));
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
	 * @param weight the weight of the edge
	 * @throws IllegalArgumentException when the specified weight is a negative number
	 */
	public void setWeight(double weight) throws IllegalArgumentException
	{
		if (checkWeight(weight))
		{
			this.weight = weight;
		}
	}
	
	/**
	 * Internal method to check whether the weight value is appropriate.
	 * @throws IllegalArgumentException when the specified weight is a negative number
	 * @return whether or not the weight is valid
	 */
	private boolean checkWeight(double weight) throws IllegalArgumentException
	{
		// TODO v2.1 log this event and continue instead of throwing an exception?
		if (weight < 0.0)
		{
			String errormsg = "The edge weight should be positive!";
			throw new IllegalArgumentException(errormsg);
		}
		return true;
	}
	
	
	@Override
	public String toString()
	{
		return EdgeIO.writeDefinitionToTab(this);
	}


}
