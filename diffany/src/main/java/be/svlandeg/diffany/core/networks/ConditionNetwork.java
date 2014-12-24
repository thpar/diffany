package be.svlandeg.diffany.core.networks;

import java.util.HashSet;
import java.util.Set;


/**
 * A kind of network that is condition-specific. There can be an unlimited number of conditions specified, but there should be at least one.
 * 
 * @author Sofie Van Landeghem
 */
public class ConditionNetwork extends InputNetwork
{

	protected Set<Condition> conditions;

	/**
	 * Create a new condition-specific network.
	 * 
	 * @param conditions at least 1 condition describing the experimental conditions  (not null or empty!)
	 * @param name the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network - can be left empty or null
	 * 
	 * @throws IllegalArgumentException when the conditions are null or empty
	 * 
	 */
	public ConditionNetwork(String name, int ID, Set<Attribute> nodeAttributes, Set<Condition> conditions) throws IllegalArgumentException
	{
		super(name, ID, nodeAttributes);
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
	 * @param name the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network - can be left empty or null
	 * @param condition one condition describing the experimental condition  (not null)
	 * 
	 * @throws IllegalArgumentException when the conditions are null or empty
	 */
	public ConditionNetwork(String name, int ID, Set<Attribute> nodeAttributes, Condition condition) throws IllegalArgumentException
	{
		super(name, ID, nodeAttributes);
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
