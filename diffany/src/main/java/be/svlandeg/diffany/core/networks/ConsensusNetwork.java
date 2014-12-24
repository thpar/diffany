package be.svlandeg.diffany.core.networks;

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

import java.util.Set;


/**
 * A consensus network is the counterpart of a differential network: containing everything but the differential edges 
 * between 2 (or more) networks, a consensus network stores the consensus or common edges between them.
 * Such a network may be useful to study 'house-keeping' interactions.
 * 
 * In contrast to a differential network, a consensus network can also be between two condition-specific networks 
 * and doesn't necessarily need a reference network.
 * 
 * @author Sofie Van Landeghem
 */
public class ConsensusNetwork extends Network
{
	
	protected Set<Network> originalNetworks;

	/**
	 * Create a new consensus network, referring to the original networks it was created from.
	 * 
	 * @param name the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param originalNetworks the original networks (at least 2!)
	 * 
	 * @throws IllegalArgumentException when the list of original networks is null or contains less than 2 networks
	 */
	public ConsensusNetwork(String name, int ID, Set<Network> originalNetworks) throws IllegalArgumentException
	{
		super(name, ID, null);
		if (originalNetworks == null || originalNetworks.size() < 1)
		{
			String errormsg = "Please define at least 1 original network!";
			throw new IllegalArgumentException(errormsg);
		}
		defineCommonAttributes(originalNetworks);
		this.originalNetworks = originalNetworks;
	}
	
	@Override
	public String getStringRepresentation()
	{
		String result = name + ": consensus network between ";
		for (Network n : originalNetworks)
		{
			result += n.getName() + " / ";
		}
		result = result.substring(0, result.length() - 2);
		return result;
	}
}
