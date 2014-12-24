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
import java.util.HashMap;
import java.util.Map;

import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;
import be.svlandeg.diffany.core.visualstyle.EdgeStyle.ArrowHead;


/**
 * This abstract class tells the GUI how to draw certain edges, depending on definitions from an {@link TreeEdgeOntology}.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class TreeEdgeDrawing extends EdgeDrawing
{

	protected TreeEdgeOntology teo;
	protected Map<String, Color> parentCatToColor;
	protected Map<String, ArrowHead> parentCatToArrowHead;
	
	
	/**
	 * Create a new EdgeDrawing object, given a specific TreeEdgeOntology object.
	 * @param teo the TreeEdgeOntology that can semantically interpret the interaction types through a tree structure
	 */
	public TreeEdgeDrawing(TreeEdgeOntology teo)
	{
		this.teo = teo;
		parentCatToColor = new HashMap<String, Color>();
		parentCatToArrowHead = new HashMap<String, ArrowHead>();
	}
	
	/**
	 * Assign a specific paint object to a source category (and its children).
	 * 
	 * @param parentCat a category (also representing its children)
	 * @param p the Color object specifying its visual properties
	 * @param source whether the category is supposed to be a source category
	 * @throws IllegalArgumentException when the either of the arguments are null, when the type was previously assigned to a paint object, 
	 * or when the (source) type is not defined in this ontology
	 */
	public void addColor(String parentCat, Color p, boolean source) throws IllegalArgumentException
	{
		if (parentCat == null || p == null)
		{
			String errormsg = "The provided parent category or the paint object should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		if (parentCatToColor.containsKey(parentCat))
		{
			String errormsg = "The provided parent category ('" + parentCat + "') already has a mapped paint object!";
			throw new IllegalArgumentException(errormsg);
		}
		if (source && !teo.isDefinedSourceCat(parentCat))
		{
			String errormsg = "The provided parent category ('" + parentCat + "') is not defined in the edge ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		parentCatToColor.put(parentCat, p);
	}

	/**
	 * Assign a specific arrowhead to a source category (and its children).
	 * 
	 * @param parentCat a category (also representing its children)
	 * @param p the ArrowHead object specifying its visual properties
	 * @param source whether the category is supposed to be a source category
	 * @throws IllegalArgumentException when the either of the arguments are null, when the type was previously assigned to a paint object, or when the type is not defined in this ontology
	 */
	public void addArrowHead(String parentCat, ArrowHead p, boolean source) throws IllegalArgumentException
	{
		if (parentCat == null || p == null)
		{
			String errormsg = "The provided parent category or the ArrowHead object should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		if (parentCatToArrowHead.containsKey(parentCat))
		{
			String errormsg = "The provided parent category ('" + parentCat + "') already has a mapped ArrowHead object!";
			throw new IllegalArgumentException(errormsg);
		}
		if (source && !teo.isDefinedSourceCat(parentCat))
		{
			String errormsg = "The provided parent category ('" + parentCat + "') is not defined in the edge ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		parentCatToArrowHead.put(parentCat, p);
	}

}
