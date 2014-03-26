package be.svlandeg.diffany.networks;

import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * A kind of network that is condition-specific. There can be an unlimited number of conditions specified, but there should be at least one.
 * 
 * @author Sofie Van Landeghem
 */
public class ConditionNetwork extends Network
{

	protected Set<Condition> conditions;

	/**
	 * Create a new condition-specific network.
	 * 
	 * @param conditions at least 1 condition describing the experimental conditions  (not null or empty!)
	 * @param name the name of this network
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 * 
	 * @throws IllegalArgumentException when the conditions are null or empty
	 * 
	 */
	public ConditionNetwork(String name, Set<Condition> conditions, NodeMapper nm) throws IllegalArgumentException
	{
		super(name, nm);
		if (conditions == null || conditions.isEmpty())
		{
			String errormsg = "Please define at least 1 condition!";
			throw new IllegalArgumentException(errormsg);
		}
		this.conditions = conditions;
	}

	/**
	 * Create a new condition-specific network.
	 * 
	 * @param condition one condition describing the experimental condition  (not null)
	 * @param name the name of this network
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 * 
	 * @throws IllegalArgumentException when the conditions are null or empty
	 */
	public ConditionNetwork(String name, Condition condition, NodeMapper nm) throws IllegalArgumentException
	{
		super(name, nm);
		if (condition == null)
		{
			String errormsg = "Please define a non-null condition!";
			throw new IllegalArgumentException(errormsg);
		}
		this.conditions = new HashSet<Condition>();
		conditions.add(condition);
	}

	/**
	 * Get the set of conditions describing the experimental environment of this network.
	 * 
	 * @return the set of conditions (should not be null or empty)
	 */
	public Set<Condition> getConditions()
	{
		return conditions;
	}

	@Override
	public String getStringRepresentation()
	{
		// TODO: check length of resulting string, define cut-off?
		String result = name + ": ";
		if (conditions != null && !conditions.isEmpty())
		{
			for (Condition c : conditions)
			{
				result += c.getDescription() + " / ";
			}
			result = result.substring(0, result.length() - 2);
		}
		else
		{
			result += "no conditions defined";
		}

		return result;
	}

}
