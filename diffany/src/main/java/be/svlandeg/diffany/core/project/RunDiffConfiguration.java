package be.svlandeg.diffany.core.project;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;


/**
 * A RunDiffConfiguration defines the necessary networks needed as input for the differential algorithms specifically,
 * and should always be used in the context of a bigger {@link Project}.
 * 
 * It should contain exactly 1 reference network, at least 1 condition-specific network, 
 * and it may contain 1 or more differential network pairs.
 * 
 * A runconfiguration object should not change after construction!
 * 
 * @author Sofie Van Landeghem
 */
public class RunDiffConfiguration extends RunConfiguration
{
	
	protected ReferenceNetwork reference;
	protected Set<ConditionNetwork> conditions;
	
	
	/**
	 * Create a new configuration with a reference network and a set of condition-specific networks. 
	 * The output result set is initialized to be empty.
	 * 
	 * @param reference the reference network (not null!)
	 * @param conditions the condition-specific networks (at least 1!)
	 * 
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public RunDiffConfiguration(ReferenceNetwork reference, Set<ConditionNetwork> conditions)
	{
		super(defineInputNetworks(reference, conditions));
		setReference(reference);
		setConditions(conditions);
	}
	
	/**
	 * Create a new configuration with a reference network and a set of condition-specific networks. 
	 * The output result set is initialized to be empty.
	 * 
	 * @param reference the reference network (not null!)
	 * @param conditions the condition-specific networks (at least 1!)
	 * @param overlapNo_cutoff the number of input networks that need to overlap to be included in the overlapping network. 
	 * It will be enforced that the reference network is always one of the supporting networks.
	 * 
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public RunDiffConfiguration(ReferenceNetwork reference, Set<ConditionNetwork> conditions, int overlapNo_cutoff)
	{
		super(defineInputNetworks(reference, conditions), overlapNo_cutoff, true);
		setReference(reference);
		setConditions(conditions);
	}

	/**
	 * Join the reference network and the condition-specific networks into one set of generic networks.
	 * 
	 * param reference the reference network (not null!)
	 * @param conditions the condition-specific networks (at least 1!)
	 * @return a merged set of generic networks
	 */
	private static Set<InputNetwork> defineInputNetworks(ReferenceNetwork reference, Set<ConditionNetwork> conditions)
	{
		Set<InputNetwork> allNetworks = new HashSet<InputNetwork>();
		allNetworks.add(reference);
		allNetworks.addAll(conditions);
		return allNetworks;
	}
	
	/**
	 * Set the condition-specific networks in this configuration.
	 * By design, this is not meant to be a public method because the configuration should not change after construction.
	 * 
	 * @param conditions the condition-specific networks (not null or empty!)
	 * @throws IllegalArgumentException if the set is null or empty
	 */
	private void setConditions(Set<ConditionNetwork> conditions) throws IllegalArgumentException
	{
		if (conditions == null || conditions.isEmpty())
		{
			String errormsg = "The set of condition-specific networks should not be null or empty!";
			throw new IllegalArgumentException(errormsg);
		}
		this.conditions = conditions;
	}
	
	/**
	 * Set the reference network of this configuration.
	 * By design, this is not meant to be a public method because the configuration should not change after construction.
	 * 
	 * @param reference the reference network
	 * @throws IllegalArgumentException if the network is null
	 */
	private void setReference(ReferenceNetwork reference) throws IllegalArgumentException
	{
		if (reference == null)
		{
			String errormsg = "The reference network should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.reference = reference;
	}
	
	/**
	 * Get the reference network of this configuration, against which the condition dependent network(s) will be compared to.
	 * 
	 * @return the reference network in this configuration (should not be null)
	 */
	public ReferenceNetwork getReferenceNetwork()
	{
		return reference;
	}

	/**
	 * Get the condition-dependent network(s): 1 or many. 
	 * Should there be 0 condition-dependent networks, the complete configuration is invalid/incomplete.
	 * 
	 * @return the condition-dependent networks in this configuration (1 or more, never null or empty))
	 */
	public Collection<ConditionNetwork> getConditionNetworks()
	{
		return conditions;
	}

}
