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

import java.util.HashSet;
import java.util.Set;


/**
 * A kind of network that is condition-specific. There can be an unlimited number of conditions specified, but there should be at least one.
 * 
 * @author Sofie Van Landeghem
 */
public class ConditionNetwork extends InputNetwork
{

	protected Set<Condition> conditions;

	/**
	 * Create a new condition-specific network.
	 * 
	 * @param conditions at least 1 condition describing the experimental conditions  (not null or empty!)
	 * @param name the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network - can be left empty or null
	 * 
	 * @throws IllegalArgumentException when the conditions are null or empty
	 * 
	 */
	public ConditionNetwork(String name, int ID, Set<Attribute> nodeAttributes, Set<Condition> conditions) throws IllegalArgumentException
	{
		super(name, ID, nodeAttributes);
		if (conditions == null || conditions.isEmpty())
		{
			String errormsg = "Please define at least 1 condition!";
			throw new IllegalArgumentException(errormsg);
		}
		this.conditions = conditions;
	}

	/**
	 * Create a new condition-specific network.
	 * 
	 * @param name the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network - can be left empty or null
	 * @param condition one condition describing the experimental condition  (not null)
	 * 
	 * @throws IllegalArgumentException when the conditions are null or empty
	 */
	public ConditionNetwork(String name, int ID, Set<Attribute> nodeAttributes, Condition condition) throws IllegalArgumentException
	{
		super(name, ID, nodeAttributes);
		if (condition == null)
		{
			String errormsg = "Please define a non-null condition!";
			throw new IllegalArgumentException(errormsg);
		}
		this.conditions = new HashSet<Condition>();
		conditions.add(condition);
	}

	/**
	 * Get the set of conditions describing the experimental environment of this network.
	 * 
	 * @return the set of conditions (should not be null or empty)
	 */
	public Set<Condition> getConditions()
	{
		return conditions;
	}

	@Override
	public String getStringRepresentation()
	{
		String result = name + ": ";
		if (conditions != null && !conditions.isEmpty())
		{
			for (Condition c : conditions)
			{
				result += c.getDescription() + " / ";
			}
			result = result.substring(0, result.length() - 2);
		}
		else
		{
			result += "no conditions defined";
		}

		return result;
	}

}
