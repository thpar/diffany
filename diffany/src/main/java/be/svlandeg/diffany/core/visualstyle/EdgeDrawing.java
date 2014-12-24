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

import be.svlandeg.diffany.core.semantics.EdgeOntology;
import be.svlandeg.diffany.core.visualstyle.EdgeStyle.ArrowHead;

/**
 * This abstract class tells the GUI how to draw certain edges, depending on definitions from an {@link EdgeOntology}.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class EdgeDrawing
{
	
	/**
	 * Create a new EdgeDrawing object
	 */
	public EdgeDrawing()
	{
	}

	/**
	 * Define the Color object of an edge by edge category.
	 * 
	 * @param edgeType the type of the edge interaction
	 * 
	 * @return the color of the edge
	 * @throws IllegalArgumentException when the provided category is null or undefined in the edgeOntology
	 */
	protected abstract Color getEdgeColor(String edgeType) throws IllegalArgumentException;
	
	/**
	 * Define the ArrowHead object of an edge by edge category.
	 * 
	 * @param edgeType the type of the edge interaction
	 * 
	 * @return the arrowhead of the edge
	 * @throws IllegalArgumentException when the provided category is null or undefined in the edgeOntology
	 */
	protected abstract ArrowHead getEdgeArrowHead(String edgeType) throws IllegalArgumentException;
	
	
	/**
	 * Define the full visual style of an edge by edge category.
	 * 
	 * @param edgeType the type of the edge interaction
	 * 
	 * @return a EdgeStyle object which specifies how the edge should be drawn
	 * @throws IllegalArgumentException when the provided category is null or undefined in the edgeOntology
	 */
	public EdgeStyle getEdgeStyle(String edgeType) throws IllegalArgumentException
	{
		return new EdgeStyle(getEdgeColor(edgeType), getEdgeArrowHead(edgeType));
	}
	
	/** 
	 * For visualisation purposes, weights need to be scaled up to a certain maximum.
	 * Beyond that maximum, all edge weights (ie edge thickness) will be the same
	 *  
	 * @return the upper boundary of weights used to define the maximal edge thickness
	 */
	public abstract double getMaxWeight();
	
	/** 
	 * For visualisation purposes, weights need to be scaled from a certain minimum.
	 * From that minimum on (INCLUSIVE), all edges will be considered to be non-existing or void.
	 *  
	 * @return the lower boundary of weights used to define the minimum edge thickness (ie. 0)
	 */
	public abstract double getMinWeight();

}
