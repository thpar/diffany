package be.svlandeg.diffany.concepts;

import java.util.Set;

/**
 * A differential network only contains differential edges between 2 (or more) networks,
 * one of which is always a static reference network.
 * 
 * @author Sofie Van Landeghem
 *
 */
public class DifferentialNetwork extends Network
{
	
	protected ReferenceNetwork reference;
	protected Set<ConditionNetwork> conditionNetworks;
	
	/**
	 * Create a new differential network, referring to the static reference network
	 * and 1 or more condition-specific networks
	 * 
	 * @param name the name of this network
	 * @param reference the corresponding reference network
	 * @param conditionNetworks the corresponding condition-specific networks
	 * @throws IllegalArgumentException when the list of condition-specific networks is null or empty,
	 * or when the reference network is null
	 */
	public DifferentialNetwork(String name, ReferenceNetwork reference, Set<ConditionNetwork> conditionNetworks) 
			throws IllegalArgumentException
	{
		super(name);
		if (reference == null)
		{
			String errormsg = "Please define at least 1 reference network!";
			throw new IllegalArgumentException(errormsg);
		}
		this.reference = reference;
		if (conditionNetworks == null || conditionNetworks.isEmpty())
		{
			String errormsg = "Please define at least 1 condition-specific network!";
			throw new IllegalArgumentException(errormsg);
		}
		this.conditionNetworks = conditionNetworks;
	}

	/* (non-Javadoc)
	 * @see be.svlandeg.diffany.concepts.Network#getStringRepresentation()
	 */
	@Override
	public String getStringRepresentation()
	{
		String result = name + ": differential network between ";
		result += reference.getName();
		result += " and ";
		for (ConditionNetwork c : conditionNetworks)
		{
			result += c.getName() + " / ";
		}
		result = result.substring(0, result.length() - 2);
		return result;
	}

}
