package be.svlandeg.diffany.core.networks.merged;

import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.core.networks.Condition;
import be.svlandeg.diffany.core.networks.EdgeDefinition;

public class MergedEdgeDefinition extends EdgeDefinition
{
	
	protected Set<Condition> conditions;
	protected boolean inReference;
	protected int support;
	
	/**
	 * Create a new 'merged' edge definition from a 'normal' edge definition
	 * 
	 * @param def the definition that contains the type, symmetry/negation status, and the weight
	 * @param conditions at least 1 condition describing the experimental conditions  (not null or empty!)
	 * @param support the number of supporting networks for this edge
	 * @param inReference whether or not this edge is present in the reference network
	 * 
	 * @throws IllegalArgumentException when the specified weight is a negative number
	 */
	public MergedEdgeDefinition(EdgeDefinition def, Set<Condition> conditions, int support, boolean inReference) throws IllegalArgumentException
	{
		this(def.getType(), def.isSymmetrical(), def.getWeight(), def.isNegated(), conditions, support, inReference);
	}
	
	/**
	 * Create a new 'merged' edge definition with a certain type, direction, weight and directionality
	 * 
	 * @param type the interaction type of this edge
	 * @param symmetrical defines whether the edge is symmetrical or directed
	 * @param weight the weight or confidence of this edge (should be positive)
	 * @param negated defines whether or not the edge is negated (e.g. does NOT bind)
	 * @param conditions at least 1 condition describing the experimental conditions  (not null or empty!)
	 * @param support the number of supporting networks for this edge
	 * @param inReference whether or not this edge is present in the reference network
	 * 
	 * @throws IllegalArgumentException when the specified weight is a negative number
	 */
	public MergedEdgeDefinition(String type, boolean symmetrical, double weight, boolean negated, Set<Condition> conditions, int support, boolean inReference) throws IllegalArgumentException
	{
		super(type, symmetrical, weight, negated);
		setConditions(conditions, support);
		this.inReference = inReference;
	}

	/**
	 * Cloning constructor
	 * 
	 * @param old the MergedEdgeDefinition to be cloned
	 */
	public MergedEdgeDefinition(MergedEdgeDefinition old)
    {
	    super(old);
	    this.conditions = new HashSet<Condition>();
	    for (Condition c : old.conditions)
	    {
	    	this.conditions.add(new Condition(c));
	    }
	    this.support = old.support;
		this.inReference = old.inReference;
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
	 * TODO javadoc
	 * @return conditions
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
		result += " (" + getSupport();
		if (inReferenceNetwork())
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
