package be.svlandeg.diffany.concepts;

import java.text.DecimalFormat;

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
	
	protected static String VOID_TYPE = "*notype*";
	protected static String DEFAULT_TYPE = "unspecified_connection";
	protected static double DEFAULT_WEIGHT = 1.0;
	protected static boolean DEFAULT_SYMM = true;
	protected static boolean DEFAULT_NEG = false;
	
	protected String type;
	
	protected boolean symmetrical;
	
	protected double weight;	
	protected boolean negated;
	
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
	 * Create a new edge definition. It will be initalized to default values 
	 * ("unspecified_connection", weight 1, symmetrical, not negated).
	 */
	public EdgeDefinition()
	{
		this(DEFAULT_TYPE , DEFAULT_SYMM, DEFAULT_WEIGHT, DEFAULT_NEG);
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
	 * Obtain a void edge, for the purpose of being able to compare it to existing edges.
	 * @return a void edge (weight == 0, symmetrical == true)
	 */
	public static EdgeDefinition getVoidEdge()
	{
		return new EdgeDefinition(VOID_TYPE , DEFAULT_SYMM, 0, DEFAULT_NEG);
	}
	
	/**
	 * Obtain a default edge, equal to the one assigned when calling an empty constructor.
	 * @return a void edge (weight == 1, symmetrical == true)
	 */
	public static EdgeDefinition getDefaultEdge()
	{
		return new EdgeDefinition(DEFAULT_TYPE , DEFAULT_SYMM, DEFAULT_WEIGHT, DEFAULT_NEG);
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
		// TODO properly log this event!
		if (checkWeight(weight))
		{
			this.weight = weight;
		}
		weight = DEFAULT_WEIGHT;
	}
	
	/**
	 * Internal method to check whether the weight value is appropriate.
	 * @throws IllegalArgumentException when the specified weight is a negative number
	 * @return whether or not the weight is valid
	 */
	private boolean checkWeight(double weight) throws IllegalArgumentException
	{
		if (weight < 0.0)
		{
			//String errormsg = "The edge weight should be positive!";
			//throw new IllegalArgumentException(errormsg);
			return false;
		}
		return true;
	}
	
	/**
	 * Get a string representation of this edge definition.
	 * More specifically, print it as: edge.type - symmetrical - weight - negated.
	 * @return a string representation of this edge, ready for printing
	 */
	public String writeToTab()
	{
		DecimalFormat df = new DecimalFormat("#.##");
		String symm = "symmetrical";
		if (! symmetrical)
		{
			symm = "directed";
		}
		String neg = "negated";
		if (! negated)
		{
			neg = "not negated";
		}
		String result = type + '\t' + symm + '\t' + df.format(weight) + '\t' + neg;
		
		return result;
	}
	
	@Override
	public String toString()
	{
		return writeToTab();
	}


}
