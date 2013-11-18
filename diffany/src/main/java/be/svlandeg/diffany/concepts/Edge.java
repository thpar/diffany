package be.svlandeg.diffany.concepts;

/**
 * Abstract class that represents an edge in a network: an edge has a source and target node
 * and can have a certain weight.
 * 
 * @author Sofie Van Landeghem
 */
public class Edge
{
	
	protected String type;
	
	protected Node source;
	protected Node target;
	protected boolean symmetrical;
	
	protected double weight;	
	protected boolean negated;
	
	/**
	 * Create a new edge with specified source and target nodes, direction and weight
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 * @param weight the weight or confidence of this edge (should be positive)
	 * @param negated defines whether or not the edge is negated (e.g. does NOT bind)
	 * @throws IllegalArgumentException when the specified weight is a negative number
	 */
	public Edge(String type, Node source, Node target, boolean symmetrical, double weight, boolean negated) throws IllegalArgumentException
	{
		checkWeight(weight);
		this.type = type;
		this.source = source;
		this.target = target;
		this.symmetrical = symmetrical;
		this.weight = weight;
	}
	
	/**
	 * Create a new edge with default weight of 1.0
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 * @param negated defines whether or not the edge is negated (e.g. does NOT bind)
	 */
	public Edge(String type, Node source, Node target, boolean symmetrical, boolean negated)
	{
		this(type, source, target, symmetrical, 1.0, negated);
	}
	
	/**
	 * Create a new edge with default weight of 1.0 and negation off
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 */
	public Edge(String type, Node source, Node target, boolean symmetrical)
	{
		this(type, source, target, symmetrical, 1.0, false);
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
	 * Get the source node of this edge
	 * When the edge is symmetrical, source and target will be defined consistently, 
	 * though they can be used interchangeably in upstream code.
	 * @return the source node of this edge
	 */
	public Node getSource()
	{
		return source;
	}
	
	/**
	 * Get the target node of this edge
	 * When the edge is symmetrical, source and target will be defined consistently, 
	 * though they can be used interchangeably in upstream code.
	 * @return the target node of this edge
	 */
	public Node getTarget()
	{
		return target;
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
	
	/**
	 * Get a string representation of this edge.
	 * More specifically, print it as: source.name - target.name - edge.type - symmetrical - weight - negated.
	 * @return a string representation of this edge, ready for printing
	 */
	public String writeToTab()
	{
		String result = source.getName() + '\t' + target.getName() + '\t' + type + '\t' + symmetrical + '\t' + weight + '\t' + negated;
		return result;
	}

}
