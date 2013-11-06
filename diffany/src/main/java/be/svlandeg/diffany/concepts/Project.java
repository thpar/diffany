package be.svlandeg.diffany.concepts;

import java.util.Collection;

/**
 * A project consists of a selection of networks and all analyses performed on these networks
 * within the current session.
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
	 * @return the reference network in this project
	 */
	public abstract Network getReferenceNetwork();
	
	/**
	 * Get the condition-dependent network(s): 1 or many. 
	 * Should there be 0 condition-dependent networks, the complete project is invalid/incomplete.
	 * 
	 * @return the condition-dependent networks in this project (1 or more)
	 */
	public abstract Collection<Network> getConditionNetworks();
	
	
	
	

}
