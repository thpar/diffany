package be.svlandeg.diffany.concepts;

import java.awt.Paint;

/**
 * This class defines the visual properties of an edge, 
 * except for the edge thickness which is directly derived from the edge weight.
 * 
 * @author Sofie Van Landeghem
 */
public class VisualEdgeStyle
{
	
	private Paint paint;
	
	/**
	 * Constructor: creates a visual edge style which cannot be modified
	 * @param paint the paint object (e.g. color)
	 */
	public VisualEdgeStyle(Paint paint)
	{
		this.paint = paint;
	}
	
	/**
	 * Retrieve the paint object
	 * @return the paint object of this visual style
	 */
	public Paint getPaint()
	{
		return paint;
	}

	
}
