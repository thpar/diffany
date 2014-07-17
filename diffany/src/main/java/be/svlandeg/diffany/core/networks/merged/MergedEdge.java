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

	protected Set<Condition> conditions;
	protected boolean inReference;
	protected int support;

	/**
	 * Create a new node from a certain definition and specifying source and target nodes.
	 * The EdgeDefinition object is not kept as such, its fields are copied (to make sure there is no dependency)
	 * 
	 * @param source the source node
	 * @param target the target node
	 * @param def the edge definition specifying the type, weight, symmetry and negation of the edge
	 * @param conditions at least 1 condition describing the experimental conditions  (not null or empty!)
	 * @param inReference whether or not this edge is present in the reference network
	 */
	public MergedEdge(Node source, Node target, EdgeDefinition def, Set<Condition> conditions, int support, boolean inReference)
	{
		super(source, target, def);
		setConditions(conditions, support);
		this.inReference = inReference;
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
	 * @param inReference whether or not this edge is present in the reference network
	 * @throws IllegalArgumentException when the specified weight is a negative number
	 */
	public MergedEdge(String type, Node source, Node target, boolean symmetrical, double weight, boolean negated, Set<Condition> conditions, int support, boolean inReference) throws IllegalArgumentException
	{
		super(type, source, target, symmetrical, weight, negated);
		setConditions(conditions, support);
		this.inReference = inReference;
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
	 */
	public MergedEdge(String type, Node source, Node target, boolean symmetrical, boolean negated, Set<Condition> conditions, int support)
	{
		super(type, source, target, symmetrical, negated);
		setConditions(conditions, support);
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
	 */
	public MergedEdge(String type, Node source, Node target, boolean symmetrical, double weight, Set<Condition> conditions, int support)
	{
		super(type, source, target, symmetrical, weight);
		setConditions(conditions, support);
	}

	/**
	 * Create a new edge with default weight of 1.0 and negation off
	 * 
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 * @param conditions at least 1 condition describing the experimental conditions  (not null or empty!)
	 */
	public MergedEdge(String type, Node source, Node target, boolean symmetrical, Set<Condition> conditions, int support)
	{
		super(type, source, target, symmetrical);
		setConditions(conditions, support);
	}

	/**
	 * Define the set of conditions in which this edge is present
	 * @param conditions at least 1 condition describing the experimental conditions  (not null or empty!)
	 */
	private void setConditions(Set<Condition> conditions, int support)
	{
		if (conditions == null || conditions.isEmpty())
		{
			String errormsg = "Please define at least 1 condition!";
			throw new IllegalArgumentException(errormsg);
		}
		this.conditions = conditions;
		if (support < 0)
		{
			String errormsg = "The support needs to be a positive integer!";
			throw new IllegalArgumentException(errormsg);
		}
		this.support = support;
	}

	/**
	 * 
	 * @return
	 */
	public Set<Condition> getConditions()
	{
		return conditions;
	}

	public boolean inReferenceNetwork()
	{
		return inReference;
	}

	public int getSupport()
	{
		return support;
	}

	@Override
	public String toString()
	{
		String result = super.toString();
		result += " (" + support;
		if (inReference)
		{
			result += ", including reference network)";
		}
		else
		{
			result += ", but not in reference network)";
		}
		return result;
	}

}
