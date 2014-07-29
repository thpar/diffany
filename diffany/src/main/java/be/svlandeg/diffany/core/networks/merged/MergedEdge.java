package be.svlandeg.diffany.core.networks.merged;

import java.util.Set;

import be.svlandeg.diffany.core.networks.Condition;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.EdgeDefinition;
import be.svlandeg.diffany.core.networks.Node;

/**
 * Class that represents an edge in a {@link MergedInputNetwork}: 
 * on top of having normal {@link Edge} properties, this edge also keeps track of the conditions in which it is present,
 * the number of (input) networks that it supports, and whether or not it is present in the reference network.
 * 
 * TODO javadoc
 * 
 * @author Sofie Van Landeghem
 */
public class MergedEdge extends Edge
{

	/**
	 * Create a new node from a certain definition and specifying source and target nodes.
	 * The EdgeDefinition object is not kept as such, its fields are copied (to make sure there is no dependency)
	 * 
	 * @param source the source node
	 * @param target the target node
	 * @param mergedDef the edge definition specifying the type, weight, symmetry and negation of the edge, as well as the type of support
	 */
	public MergedEdge(Node source, Node target, MergedEdgeDefinition mergedDef)
	{
		super(source, target, mergedDef);
	}

	/**
	 * Create a new edge with specified source and target nodes, direction and weight
	 * 
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 * @param weight the weight or confidence of this edge (should be positive)
	 * @param negated defines whether or not the edge is negated (e.g. does NOT bind)
	 * @param conditions at least 1 condition describing the experimental conditions  (not null or empty!)
	 * @param support the number of supporting networks for this edge
	 * @param inReference whether or not this edge is present in the reference network
	 * @throws IllegalArgumentException when the specified weight is a negative number
	 */
	public MergedEdge(String type, Node source, Node target, boolean symmetrical, double weight, boolean negated, Set<Condition> conditions, int support, boolean inReference) throws IllegalArgumentException
	{
		this(source, target, new MergedEdgeDefinition(type, symmetrical, weight, negated, conditions, support, inReference));
	}

	/**
	 * Create a new edge with default weight of 1.0
	 * 
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 * @param negated defines whether or not the edge is negated (e.g. does NOT bind)
	 * @param conditions at least 1 condition describing the experimental conditions  (not null or empty!)
	 * @param support the number of supporting networks for this edge
	 * @param inReference whether or not this edge is present in the reference network
	 */
	public MergedEdge(String type, Node source, Node target, boolean symmetrical, boolean negated, Set<Condition> conditions, int support, boolean inReference)
	{
		this(source, target, new MergedEdgeDefinition(type, symmetrical, EdgeDefinition.DEFAULT_WEIGHT, negated, conditions, support, inReference));
	}

	/**
	 * Create a new edge, which will be defined as not negated.
	 * 
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 * @param weight the weight or confidence of this edge (should be positive)
	 * @param conditions at least 1 condition describing the experimental conditions  (not null or empty!)
	 * @param support the number of supporting networks for this edge
	 * @param inReference whether or not this edge is present in the reference network
	 */
	public MergedEdge(String type, Node source, Node target, boolean symmetrical, double weight, Set<Condition> conditions, int support, boolean inReference)
	{
		this(source, target, new MergedEdgeDefinition(type, symmetrical, weight, EdgeDefinition.DEFAULT_NEG, conditions, support, inReference));
	}

	/**
	 * Create a new edge with default weight of 1.0 and negation off
	 * 
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 * @param conditions at least 1 condition describing the experimental conditions  (not null or empty!)
	 * @param support the number of supporting networks for this edge
	 * @param inReference whether or not this edge is present in the reference network
	 */
	public MergedEdge(String type, Node source, Node target, boolean symmetrical, Set<Condition> conditions, int support, boolean inReference)
	{
		this(source, target, new MergedEdgeDefinition(type, symmetrical, EdgeDefinition.DEFAULT_WEIGHT, EdgeDefinition.DEFAULT_NEG, conditions, support, inReference));
	}

	/**
	 * TODO javadoc
	 * @return conditions
	 */
	public Set<Condition> getConditions()
	{
		return ((MergedEdgeDefinition) def).conditions;
	}

	public boolean inReferenceNetwork()
	{
		return ((MergedEdgeDefinition) def).inReference;
	}

	public int getSupport()
	{
		return ((MergedEdgeDefinition) def).support;
	}

}
