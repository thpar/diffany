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


/**
 * This class generates special edges to be used in the input/output networks, such as default edges,
 * virtual edges or void edges.
 * 
 * @author Sofie Van Landeghem
 */
public class EdgeGenerator
{
	
	protected static String GENERIC_DIRECT_TYPE = "*generic_directed_connection*";
	protected static String GENERIC_SYMM_TYPE = "*generic_symmetrical_connection*";
	
	
	/**
	 * Obtain a default edge. It will be initalized to default values from {@link EdgeDefinition}: "unspecified_connection", weight 1, symmetrical, not negated.
	 * 
	 * @return a default edge (weight == 1, symmetrical == true)
	 */
	public EdgeDefinition getDefaultEdge()
	{
		return new EdgeDefinition(EdgeDefinition.DEFAULT_TYPE, EdgeDefinition.DEFAULT_SYMM, EdgeDefinition.DEFAULT_WEIGHT, EdgeDefinition.DEFAULT_NEG);
	}
	
	
	/**
	 * Obtain a void edge, for the purpose of being able to compare it to existing edges.
	 * 
	 * @param symmetrical whether or not the edge should be symmetrical
	 * @return a void edge (weight == 0)
	 */
	public EdgeDefinition getVoidEdge(boolean symmetrical)
	{
		if (symmetrical)
		{
			return new EdgeDefinition(GENERIC_SYMM_TYPE , true, 0, EdgeDefinition.DEFAULT_NEG);
		}
		return new EdgeDefinition(GENERIC_DIRECT_TYPE , false, 0, EdgeDefinition.DEFAULT_NEG);
	}
	
	
	

}
