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
	//private enum arrowHead{};
	
	/**
	 * Constructor: creates a visual edge style which cannot be modified
	 * @param color object 
	 */
	public VisualEdgeStyle(Color color)
	{
		this.color = color;
	}
	
	/**
	 * Retrieve the paint object
	 * @return the paint object of this visual style
	 */
	public Color getColor()
	{
		return color;
	}

	
}
