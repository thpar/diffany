package be.svlandeg.diffany.core.project;

/*
 * #%L
 * Diffany
 * %%
 * Copyright (C) 2014 PSB/UGent - Sofie Van Landeghem and Thomas Van Parys
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU General Lesser Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/lgpl-3.0.html>.
 * #L%
 */


import java.util.Collection;
import java.util.Set;

import be.svlandeg.diffany.core.networks.InputNetwork;


/**
 * A RunConfiguration defines the necessary networks needed as input for the Diffany algorithms,
 * and should always be used in the context of a bigger {@link Project}.
 * 
 * This configuration can be used to calculate consensus networks, but differential networks can only be calculated with a RunDiffConfiguration object.
 * 
 * A runconfiguration object should not change after construction!
 * 
 * @author Sofie Van Landeghem
 */
public class RunConfiguration
{
	
	protected Set<InputNetwork> inputNetworks;
	protected int supportCutoff;
	protected boolean refRequired;
	
	
	/**
	 * Create a new configuration with a set of input networks and a required supportCutoff (between 2 and the size of the network set).
	 * The output result set is initialized to be empty.
	 * 
	 * @param inputNetworks the input networks (not null or empty!)
	 * @param supportCutoff the number of input networks that need to support an edge for it to be included in the consensus network
	 * @param refRequired whether or not the presence of the edge in the reference network is required for it to be included in the consensus network
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public RunConfiguration(Set<InputNetwork> inputNetworks, int supportCutoff, boolean refRequired)
	{
		setInputs(inputNetworks);
		setSupportCutoff(supportCutoff);
		setRefRequired(refRequired);
	}
	
	/**
	 * Create a new configuration with a set of input networks. The required support cutoff is by default set to the size of this set. 
	 * The output result set is initialized to be empty.
	 * 
	 * @param inputNetworks the input networks (not null or empty!)
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public RunConfiguration(Set<InputNetwork> inputNetworks)
	{
		this(inputNetworks, inputNetworks.size(), false);
	}
	
	/**
	 * Set the input networks in this configuration. By design, this is not meant to be a public method because the configuration should not change after construction.
	 * 
	 * @param inputNetworks the input networks (not null or empty!)
	 * @throws IllegalArgumentException if the set is null or empty
	 */
	protected void setInputs(Set<InputNetwork> inputNetworks) throws IllegalArgumentException
	{
		if (inputNetworks == null || inputNetworks.isEmpty())
		{
			String errormsg = "The set of input networks should not be null or empty!";
			throw new IllegalArgumentException(errormsg);
		}
		this.inputNetworks = inputNetworks;
	}
	
	/**
	 * Alter the fact whether or not the reference network needs to have an edge for it to be allowed inclusion in the consensus network.
	 * By design, this is not meant to be a public method because the configuration should not change after construction.
	 * 
	 * @param refRequired the new boolean state
	 */
	protected void setRefRequired(boolean refRequired)
	{
		this.refRequired = refRequired;
	}
	
	/**
	 * Retrieve whether or not the reference network needs to have an edge for it to be allowed inclusion in the consensus network
	 * @return the value of the refRequired
	 */
	public boolean getRefRequired()
	{
		return refRequired;
	}
	
	/**
	 * Alter the required number of supporting networks needed before an edge will be present in the consensus network.
	 * By design, this is not meant to be a public method because the configuration should not change after construction.
	 * 
	 * @param supportCutoff the new cutoff
	 */
	protected void setSupportCutoff(int supportCutoff)
	{
		if (supportCutoff < 2)
		{
			String errormsg = "The support cutoff should at least be 2!";
			throw new IllegalArgumentException(errormsg);
		}
		this.supportCutoff = supportCutoff;
	}
	
	/**
	 * Return the required number of supporting networks needed before the edge will be present in the consensus network.
	 * When a differential network is being calculated, this threshold is reduced by one (reference network) to get 
	 * the required representation of the condition-dependent networks
	 * 
	 * @return the minimal number of required networks for consensus edges, which is the size of the input set unless specifically defined otherwise.
	 */
	public int getSupportCutoff()
	{
		return supportCutoff;
	}
	
	/**
	 * Get all input network(s): 1 or many. 
	 * 
	 * @return the input networks in this configuration (1 or more, never null or empty))
	 */
	public Collection<InputNetwork> getInputNetworks()
	{
		return inputNetworks;
	}

}
