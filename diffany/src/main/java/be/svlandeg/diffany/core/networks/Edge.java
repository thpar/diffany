package be.svlandeg.diffany.core.networks;

import be.svlandeg.diffany.core.io.EdgeIO;

/**
 * Class that represents an edge in a network: an edge has a source and target node
 * and can have a certain weight. It is symmetrical or not, and may or may not be negated.
 * 
 * Edge data can be saved and loaded through the {@link EdgeIO} class
 * 
 * @author Sofie Van Landeghem
 */
public class Edge extends EdgeDefinition
{

	protected Node source;
	protected Node target;

	/**
	 * Create a new node from a certain definition and specifying source and target nodes.
	 * The EdgeDefinition object is not kept as such, its fields are copied (to make sure there is no dependency).
	 * 
	 * @param source the source node
	 * @param target the target node
	 * @param def the edge definition specifying the type, weight, symmetry and negation of the edge
	 */
	public Edge(Node source, Node target, EdgeDefinition def)
	{
		this(def.type, source, target, def.symmetrical, def.weight, def.negated);
	}

	/**
	 * Create a new edge with specified source and target nodes, direction and weight.
	 * 
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
	 * Create a new edge with default weight of 1.
	 * 
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 * @param negated defines whether or not the edge is negated (e.g. does NOT bind)
	 */
	public Edge(String type, Node source, Node target, boolean symmetrical, boolean negated)
	{
		this(type, source, target, symmetrical, DEFAULT_WEIGHT, negated);
	}

	/**
	 * Create a new edge, which will be defined as not negated.
	 * 
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 * @param weight the weight or confidence of this edge (should be positive)
	 */
	public Edge(String type, Node source, Node target, boolean symmetrical, double weight)
	{
		this(type, source, target, symmetrical, weight, DEFAULT_NEG);
	}

	/**
	 * Create a new edge with default weight of 1.0 and negation off.
	 * 
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 */
	public Edge(String type, Node source, Node target, boolean symmetrical)
	{
		this(type, source, target, symmetrical, DEFAULT_WEIGHT, DEFAULT_NEG);
	}

	/**
	 * Get the source node of this edge.
	 * When the edge is symmetrical, source and target will be defined consistently,
	 * though they can be used interchangeably in upstream code.
	 * 
	 * @return the source node of this edge
	 */
	public Node getSource()
	{
		return source;
	}

	/**
	 * Get the target node of this edge.
	 * When the edge is symmetrical, source and target will be defined consistently,
	 * though they can be used interchangeably in upstream code.
	 * 
	 * @return the target node of this edge
	 */
	public Node getTarget()
	{
		return target;
	}

	@Override
	public String toString()
	{
		return EdgeIO.writeToTab(this);
	}
	
	/** 
	 * Retrieve whether or not an edge is virtual. It is considered virtual if the source and/or target nodes are virtual.
	 * @return whether or not this is a virtual edge
	 */
	public boolean isVirtual()
	{
		boolean virtual = false;
		if (source.isVirtual() || target.isVirtual())
		{
			virtual = true;
		}
		return virtual;
	}

}
