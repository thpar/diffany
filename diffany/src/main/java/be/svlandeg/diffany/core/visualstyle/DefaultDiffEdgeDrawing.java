package be.svlandeg.diffany.core.visualstyle;

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
		// TODO v2.1: is this reasonable?
		return max_weight;
	}

	@Override
	public double getMinWeight()
	{
		return min_weight;
	}

}
