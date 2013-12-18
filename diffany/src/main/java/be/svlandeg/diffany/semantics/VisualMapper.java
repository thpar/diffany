package be.svlandeg.diffany.semantics;

import java.awt.Color;
import java.awt.Paint;

/** 
 * This class defines the visual style of a network, and will be tightly linked to an edge ontology
 * as the colors/shapes are dependent on the underlying semantics of the nodes/edges. 
 * 
 * @author Sofie Van Landeghem
 *
 */
public class VisualMapper
{
	
	/**
	 * Define the visual style of an edge in a differential network, by edge category.
	 * 
	 * @param eo the edge ontology which can interpret the semantics of the edge
	 * @param category the category of the edge interaction
	 * @return a Paint object which specifies how the edge should be drawn
	 */
	public Paint getDifferentialEdgeStyle(EdgeOntology eo, String category)
	{
		return Color.PINK;
	}
	
	/**
	 * Define the visual style of an edge in a 'normal' network (reference, condition-dependent or overlap).
	 * 
	 * @param eo the edge ontology which can interpret the semantics of the edge
	 * @param type the type of the edge interaction
	 * @return a Paint object which specifies how the edge should be drawn
	 */
	public Paint getStaticEdgeStyle(EdgeOntology eo, String type)
	{
		return Color.ORANGE;
	}
	
}
