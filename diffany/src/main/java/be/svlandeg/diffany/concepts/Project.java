package be.svlandeg.diffany.concepts;

import java.util.Collection;

/**
 * A project consists of a selection of networks and all analyses performed on these networks
 * within the current session.
 * It should contain exactly 1 reference network, at least 1 condition-specific network,
 * and it may contain 1 or more differential networks.
 * 
 * @author Sofie Van Landeghem
 *
 */
public abstract class Project
{
	/**
	 * Get the reference network of this project, against which the condition dependent network(s)
	 * will be compared to.
	 * 
	 * @return the reference network in this project (should not be null)
	 */
	public abstract ReferenceNetwork getReferenceNetwork();
	
	/**
	 * Get the condition-dependent network(s): 1 or many. 
	 * Should there be 0 condition-dependent networks, the complete project is invalid/incomplete.
	 * 
	 * @return the condition-dependent networks in this project (1 or more, never null or empty))
	 */
	public abstract Collection<ConditionNetwork> getConditionNetworks();
	
	/**
	 * Get the differential networks in the project: 0, 1 or more 
	 *
	 * @return the differential networks in this project (if any, otherwise empty set, but never null)
	 */
	public abstract Collection<DifferentialNetwork> getDifferentialNetworks();
	
	

}
