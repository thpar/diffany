package be.svlandeg.diffany.core.project;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;


/**
 * A RunConfiguration defines the necessary networks needed as input for the Diffany algorithms,
 * and should always be used in the context of a bigger {@link Project}.
 * 
 * It should contain exactly 1 reference network, at least 1 condition-specific network, 
 * and it may contain 1 or more differential networks.
 * 
 * @author Sofie Van Landeghem
 */
public class RunConfiguration
{
	
	protected ReferenceNetwork reference;
	protected Set<ConditionNetwork> conditions;
	protected Set<DifferentialNetwork> differentials;
	
	
	/**
	 * Create a new configuration with a reference network and a set of condition-specific networks. 
	 * The set of (output) differential networks is initialized to be empty.
	 * 
	 * @param reference the reference network (not null!)
	 * @param conditions the condition-specific networks (at least 1!)
	 * 
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public RunConfiguration(ReferenceNetwork reference, Set<ConditionNetwork> conditions)
	{
		setReference(reference);
		setConditions(conditions);
		differentials = new HashSet<DifferentialNetwork>();
	}
	
	/**
	 * Set the condition-specific networks in this configuration.
	 * (currently not a public method - changes to it would influence the differential networks (TODO v3.0))
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
	 * (currently not a public method - changes to it would influence the differential networks (TODO v3.0))
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
	 * Add a differential network to this configuration.
	 * 
	 * @param differential a new differential network
	 */
	public void addDifferential(DifferentialNetwork differential)
	{
		differentials.add(differential);
	}
	
	/**
	 * Get the reference network of this configuration, against which the condition
	 * dependent network(s) will be compared to.
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

	/**
	 * Get the differential networks in the configuration: 0, 1 or more
	 * 
	 * @return the differential networks in this configuration (if any, otherwise empty set, but never null)
	 */
	public Collection<DifferentialNetwork> getDifferentialNetworks()
	{
		return differentials;
	}

}
