package be.svlandeg.diffany.concepts;

import java.util.Set;

/**
 * A differential network only contains differential edges between 2 (or more) networks,
 * one of which is always a 'static' reference network.
 * 
 * @author Sofie Van Landeghem
 */
public class DifferentialNetwork extends Network
{
	
	protected ReferenceNetwork reference;
	protected Set<ConditionNetwork> conditionNetworks;
	protected SharedNetwork shared;
	
	
	/**
	 * Create a new differential network, referring to exactly one static reference network
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
	
	/**
	 * Set the shared ('house-keeping') network that is associated to this differential network
	 * @param shared the complementary shared network
	 */
	public void setSharedNetwork(SharedNetwork shared)
	{
		this.shared = shared;
	}
	
	/**
	 * Get the shared ('house-keeping') network associated to this differential network
	 * @return the shared network that complements this differential network (may be null)
	 */
	public SharedNetwork getSharedNetwork()
	{
		return shared;
	}
	
	/**
	 * Get the reference network associated to this differential network
	 * @return the reference network 
	 */
	public ReferenceNetwork getReferenceNetwork()
	{
		return reference;
	}
	
	/**
	 * Get the (set of) condition-specific networks associated to this differential network
	 * @return the (set of) condition-specific networks
	 */
	public Set<ConditionNetwork> getConditionNetworks()
	{
		return conditionNetworks;
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
