package be.svlandeg.diffany.concepts;


/**
 * Class that represents an edge in a network: an edge has a source and target node
 * and can have a certain weight. It is symmetrical or not, and may or may not be negated.
 * 
 * @author Sofie Van Landeghem
 */
public class Edge extends EdgeDefinition
{

	protected Node source;
	protected Node target;

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
		super(type, symmetrical, weight, negated);
		this.source = source;
		this.target = target;
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
	 * Create a new edge with default weight of 1.0
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 * @param weight the weight or confidence of this edge (should be positive)
	 */
	public Edge(String type, Node source, Node target, boolean symmetrical, double weight)
	{
		this(type, source, target, symmetrical, weight, false);
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
	 * Create a new node from a certain definition and specifying source and target nodes.
	 * @param source the source node
	 * @param target the target node
	 * @param def the edge definition specifying the type, weight, symmetry and negation of the edge
	 */
	public Edge(Node source, Node target, EdgeDefinition def)
	{
		this(def.type, source, target, def.symmetrical, def.weight, def.negated);
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
	 * Get a string representation of this edge.
	 * More specifically, print it as: source.name - target.name - edge.type - symmetrical - weight - negated.
	 * @return a string representation of this edge, ready for printing
	 */
	public String writeToTab()
	{
		String defResult = super.writeToTab();
		String result = source.getName() + '\t' + target.getName() + '\t' + defResult;
		
		return result;
	}

}
