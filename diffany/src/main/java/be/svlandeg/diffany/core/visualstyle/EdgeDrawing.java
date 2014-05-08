package be.svlandeg.diffany.core.visualstyle;

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