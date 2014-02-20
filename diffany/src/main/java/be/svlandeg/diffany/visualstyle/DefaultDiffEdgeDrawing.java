package be.svlandeg.diffany.visualstyle;

import java.awt.Color;

import be.svlandeg.diffany.semantics.TreeEdgeOntology;
import be.svlandeg.diffany.visualstyle.EdgeStyle.ArrowHead;

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
	public ArrowHead getEdgeArrowHead(String category)
	{
		if (category == null)
		{
			String errormsg = "The provided differential type ('" + category + "') should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		if (teo.isPosDirected(category))
		{
			return pos_ah;
		}
		if (teo.isNegDirected(category))
		{
			return neg_ah;
		}
		if (teo.isNegSymm(category) || teo.isPosSymm(category))
		{
			return symm_ah;
		}
		return default_ah;
	}

	@Override
	public Color getEdgeColor(String category)
	{
		if (category == null)
		{
			String errormsg = "The provided differential type ('" + category + "') should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		if (teo.isPosDirected(category) || teo.isPosSymm(category))
		{
			return pos_paint;
		}
		if (teo.isNegSymm(category) || teo.isNegDirected(category))
		{
			return neg_paint;
		}
		return default_paint;
	}

	@Override
	public double getMaxWeight()
	{
		// TODO v2.0: is this reasonable?
		return max_weight;
	}

	@Override
	public double getMinWeight()
	{
		return min_weight;
	}

}
