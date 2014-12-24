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
 * This class provides a default definition of drawing edges in a GUI for the source networks. 
 * Specifically, this covers the input networks (reference and condition-dependent) as wel as the consensus networks.
 * 
 * @author Sofie Van Landeghem
 */
public class DefaultSourceEdgeDrawing extends TreeEdgeDrawing
{

	protected static ArrowHead symm_ah = ArrowHead.NONE;
	protected static ArrowHead neutral_ah = ArrowHead.NONE;

	protected static Color neutral_paint = Color.LIGHT_GRAY;
	
	protected static double min_weight = 0;
	protected static double max_weight = 20;

	/**
	 * Create a new DefaultSourceEdgeDrawing object, given a specific EdgeOntology object.
	 * 
	 * @param teo the TreeEdgeOntology that can semantically interpret the interaction types
	 */
	public DefaultSourceEdgeDrawing(TreeEdgeOntology teo)
	{
		super(teo);
	}

	@Override
	protected Color getEdgeColor(String edgeType)
	{
		if (teo.isDefinedSourceType(edgeType))
		{
			String childCat = teo.getSourceCategory(edgeType);
			Color foundColor = parentCatToColor.get(childCat);
			
			while (foundColor == null && childCat != null)
			{
				String parentCat = teo.retrieveCatParent(childCat);
				foundColor = parentCatToColor.get(parentCat);
				
				childCat = parentCat;
			}
			if (foundColor != null)
			{
				return foundColor;
			}
			return neutral_paint;
		}
		String errormsg = "The provided source type ('" + edgeType + "') is not known in the edge ontology!";
		errormsg += " Make sure to register all source networks to the Project before running the visualisation.";
		throw new IllegalArgumentException(errormsg);
	}

	@Override
	protected ArrowHead getEdgeArrowHead(String edgeType)
	{
		if (teo.isDefinedSourceType(edgeType))
		{
			String childCat = teo.getSourceCategory(edgeType);
			ArrowHead foundArrowHead = parentCatToArrowHead.get(childCat);
			while (foundArrowHead == null && childCat != null)
			{
				String parentCat = teo.retrieveCatParent(childCat);
				foundArrowHead = parentCatToArrowHead.get(parentCat);
				childCat = parentCat;
			}
			if (foundArrowHead != null)
			{
				return foundArrowHead;
			}
			if (teo.isSymmetricalSourceType(edgeType))
			{
				return symm_ah;
			}
			return neutral_ah;
		}
		String errormsg = "The provided source type ('" + edgeType + "') is not known in the edge ontology!";
		errormsg += " Make sure to register all source networks to the Project before running the visualisation.";
		throw new IllegalArgumentException(errormsg);
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
