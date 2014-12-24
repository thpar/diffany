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

import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;
import be.svlandeg.diffany.core.visualstyle.EdgeStyle.ArrowHead;

/**
 * This class provides a default definition of drawing edges in a GUI for the differential network(s).
 * 
 * @author Sofie Van Landeghem
 */
public class DefaultDiffEdgeDrawing extends TreeEdgeDrawing
{

	protected static Color neg_paint = Color.RED;
	protected static Color pos_paint = Color.GREEN;
	protected static Color default_paint = Color.GRAY;

	protected static ArrowHead neg_ah = ArrowHead.ARROW;
	protected static ArrowHead pos_ah = ArrowHead.ARROW;
	protected static ArrowHead default_ah = ArrowHead.NONE;

	protected static ArrowHead symm_ah = ArrowHead.NONE;

	protected static double min_weight = 0;
	protected static double max_weight = 20;

	/**
	 * Create a new DefaultDiffEdgeDrawing object, given a specific EdgeOntology object.
	 * 
	 * @param teo the TreeEdgeOntology that can semantically interpret the interaction types
	 */
	public DefaultDiffEdgeDrawing(TreeEdgeOntology teo)
	{
		super(teo);
	}

	@Override
	public ArrowHead getEdgeArrowHead(String type)
	{
		if (type == null)
		{
			String errormsg = "The provided differential type ('" + type + "') should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		
		// TODO v.3.0 allow tree-like traversal of differential edge types (cf. DefaultSourceEdgeDrawing.getEdgeColor)
		ArrowHead foundArrowHead = parentCatToArrowHead.get(type);
		
		if (foundArrowHead != null)
		{
			return foundArrowHead;
		}
		
		if (teo.isPosDirected(type))
		{
			return pos_ah;
		}
		if (teo.isNegDirected(type))
		{
			return neg_ah;
		}
		if (teo.isNegSymm(type) || teo.isPosSymm(type))
		{
			return symm_ah;
		}
		return default_ah;
	}

	@Override
	public Color getEdgeColor(String type)
	{
		if (type == null)
		{
			String errormsg = "The provided differential type ('" + type + "') should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		
		// TODO v.3.0 allow tree-like traversal of differential edge types (cf. DefaultSourceEdgeDrawing.getEdgeColor)
		Color foundColor = parentCatToColor.get(type);
		
		if (foundColor != null)
		{
			return foundColor;
		}
		
		if (teo.isPosDirected(type) || teo.isPosSymm(type))
		{
			return pos_paint;
		}
		if (teo.isNegSymm(type) || teo.isNegDirected(type))
		{
			return neg_paint;
		}
		return default_paint;
	}

	@Override
	public double getMaxWeight()
	{
		// TODO v2.2: is this reasonable?
		return max_weight;
	}

	@Override
	public double getMinWeight()
	{
		return min_weight;
	}

}
