package be.svlandeg.diffany.concepts;

import java.awt.Color;

/**
 * This class defines the visual properties of an edge, 
 * except for the edge thickness which is directly derived from the edge weight.
 * 
 * @author Sofie Van Landeghem
 */
public class VisualEdgeStyle
{
	
	private Color color;
	private ArrowHead ah;
	
	public enum ArrowHead{ARROW, T, NONE, DIAMOND};
	
	
	/**
	 * Constructor: creates a visual edge style which cannot be modified
	 * @param color object 
	 */
	public VisualEdgeStyle(Color color, ArrowHead ah)
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
