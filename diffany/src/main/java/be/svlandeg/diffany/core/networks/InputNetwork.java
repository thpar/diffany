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
 * A generic input network is not yet defined as a specific Diffany network type, and can at times be used to generically read input data,
 * which still requires pre-processing.
 * 
 * @author Sofie Van Landeghem
 *
 */
public class InputNetwork extends Network
{

	/**
	 * Create a new generic input network.
	 * 
	 * @param name the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network - can be left empty or null
	 * 
	 */
	public InputNetwork(String name, int ID, Set<Attribute> nodeAttributes)
	{
		super(name, ID, nodeAttributes);
	}

	/**
	 * Create a new generic input network.
	 * 
	 * @param networkName the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network
	 * @param nodes the nodes of this network
	 * @param edges the edges of this network
	 */
	public InputNetwork(String networkName, int ID, Set<Attribute> nodeAttributes, Set<Node> nodes, Set<Edge> edges)
	{
		super(networkName, ID, nodeAttributes, nodes, edges);
	}

	@Override
	public String getStringRepresentation()
	{
		return name + ": input network";
	}
}
