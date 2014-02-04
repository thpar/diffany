package be.svlandeg.diffany.visualstyle;

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
	
	// These types correspond directly to Cytoscape types, but can also be used in other visualisation tools
	public enum ArrowHead{ARROW, T, NONE, DIAMOND};
	
	
	/**
	 * Constructor: creates a visual edge style which cannot be modified
	 * @param color object 
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
