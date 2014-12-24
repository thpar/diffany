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
 * A reference network is used as comparison against condition-dependent networks.
 * It can be seen as a 'static' network with unspecified/unknown/wild-type conditions.
 * 
 * @author Sofie Van Landeghem
 *
 */
public class ReferenceNetwork extends InputNetwork
{
	
	/**
	 * Create a new static reference network.
	 * 
	 * @param name the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network - can be left empty or null
	 * 
	 */
	public ReferenceNetwork(String name, int ID, Set<Attribute> nodeAttributes)
	{
		super(name, ID, nodeAttributes);
	}

	/* (non-Javadoc)
	 * @see be.svlandeg.diffany.core.networks.Network#getStringRepresentation()
	 */
	@Override
	public String getStringRepresentation()
	{
		return name + ": input reference network";
	}

	

	
}
