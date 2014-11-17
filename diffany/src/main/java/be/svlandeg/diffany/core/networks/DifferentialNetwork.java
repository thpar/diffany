package be.svlandeg.diffany.core.networks;

import java.util.HashSet;
import java.util.Set;


/**
 * A differential network contains differential edges between 2 (or more) networks,
 * one of which is always a 'static' reference network.
 * 
 * @author Sofie Van Landeghem
 */
public class DifferentialNetwork extends Network
{
	
	protected ReferenceNetwork reference;
	protected Set<ConditionNetwork> conditionNetworks;
	
	
	/**
	 * Create a new differential network, referring to exactly one static reference network
	 * and 1 or more condition-specific networks
	 * 
	 * @param name the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param reference the corresponding reference network (not null!)
	 * @param conditionNetworks the corresponding condition-specific networks (not null or empty!)
	 * 
	 * @throws IllegalArgumentException when the list of condition-specific networks is null or empty,
	 * or when the reference network is null
	 */
	public DifferentialNetwork(String name, int ID, 
			ReferenceNetwork reference, Set<ConditionNetwork> conditionNetworks) 
			throws IllegalArgumentException
	{
		super(name, ID, null);
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
		
		Set<Network> originalNetworks = new HashSet<Network>();
		originalNetworks.add(reference);
		originalNetworks.addAll(conditionNetworks);
		defineCommonAttributes(originalNetworks);
	}
	
	/**
	 * Create a new differential network, referring to exactly one static reference network and one condition-specific network
	 * 
	 * @param name the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param reference the corresponding reference network (not null!)
	 * @param conditionNetwork the corresponding condition-specific network (not null!)
	 * 
	 * @throws IllegalArgumentException when the list of condition-specific networks is null or empty,
	 * or when the reference network is null
	 */
	public DifferentialNetwork(String name, int ID, ReferenceNetwork reference, ConditionNetwork conditionNetwork) 
			throws IllegalArgumentException
	{
		super(name, ID, null);
		if (reference == null)
		{
			String errormsg = "Please define at least 1 reference network!";
			throw new IllegalArgumentException(errormsg);
		}
		this.reference = reference;
		if (conditionNetwork == null)
		{
			String errormsg = "Please define a non-null condition-specific network!";
			throw new IllegalArgumentException(errormsg);
		}
		conditionNetworks = new HashSet<ConditionNetwork>();
		conditionNetworks.add(conditionNetwork);
		
		Set<Network> originalNetworks = new HashSet<Network>();
		originalNetworks.add(reference);
		originalNetworks.add(conditionNetwork);
		defineCommonAttributes(originalNetworks);
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
	 * @see be.svlandeg.diffany.core.networks.Network#getStringRepresentation()
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
