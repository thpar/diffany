package be.svlandeg.diffany.concepts;

import java.util.Set;

/**
 * A kind of network that is condition-specific. There can be an unlimited
 * number of conditions specified, but there should be at least one.
 * 
 * @author Sofie Van Landeghem
 */
public class ConditionNetwork extends Network
{

	protected Set<Condition> conditions;

	/**
	 * Create a new condition-specific network.
	 * 
	 * @param conditions
	 *            at least 1 condition describing the experimental conditions.
	 * @param name the name of this network
	 * @throws IllegalArgumentException
	 *             when the conditions are null or empty
	 * 
	 */
	public ConditionNetwork(String name, Set<Condition> conditions) throws IllegalArgumentException
	{
		super(name);
		if (conditions == null || conditions.isEmpty())
		{
			String errormsg = "Please define at least 1 condition!";
			throw new IllegalArgumentException(errormsg);
		}
		this.conditions = conditions;
	}

	/**
	 * Get the set of conditions describing the experimental environment of this
	 * network.
	 * 
	 * @return the set of conditions (should not be null or empty)
	 */
	public Set<Condition> getConditions()
	{
		return conditions;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * be.svlandeg.diffany.concepts.Network#getStringRepresentation()
	 */
	@Override
	public String getStringRepresentation()
	{
		// TODO: check length of resulting string, define cut-off?
		String result = name + ": "; 
		if (conditions != null && ! conditions.isEmpty())
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
