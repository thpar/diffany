package be.svlandeg.diffany.core.visualstyle;

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

import java.awt.Color;

/**
 * This class defines the main visual properties of an edge, such as the color and type of arrowhead.
 * Once defined, the properties can not be changed (create a new object instead).
 * 
 * @author Sofie Van Landeghem
 */
public class EdgeStyle
{
	
	private final Color color;
	private final ArrowHead ah;
	
	/* These types correspond directly to Cytoscape types, but can also be used in other visualisation tools */
	public enum ArrowHead{ARROW, T, NONE, DIAMOND};
	
	
	/**
	 * Constructor: creates a visual edge style which cannot be modified
	 * @param color the color of the edge 
	 * @param ah the type of arrowhead
	 */
	public EdgeStyle(Color color, ArrowHead ah)
	{
		this.color = color;
		this.ah = ah;
	}
	
	/**
	 * Retrieve the color
	 * @return the color of this visual style
	 */
	public Color getColor()
	{
		return color;
	}
	
	/**
	 * Retrieve the arrowhead style
	 * @return the type of arrowhead of this visual style
	 */
	public ArrowHead getArrowHead()
	{
		return ah;
	}

	
}
