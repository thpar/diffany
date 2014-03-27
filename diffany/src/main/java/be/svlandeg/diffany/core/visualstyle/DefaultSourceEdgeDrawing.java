package be.svlandeg.diffany.core.visualstyle;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;

import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;
import be.svlandeg.diffany.core.visualstyle.EdgeStyle.ArrowHead;

/**
 * This class provides a default definition of drawing edges in a GUI for the source networks. Specifically, this covers the input networks (reference and condition-dependent) as wel as the overlapping networks.
 * 
 * @author Sofie Van Landeghem
 */
public class DefaultSourceEdgeDrawing extends TreeEdgeDrawing
{

	protected Map<String, Color> parentCatToColor;
	protected Map<String, ArrowHead> parentCatToArrowHead;

	protected static ArrowHead symm_ah = ArrowHead.NONE;
	protected static ArrowHead neutral_ah = ArrowHead.NONE;

	protected static Color neutral_paint = Color.LIGHT_GRAY;
	
	protected static double min_weight = 0;
	protected static double max_weight = 20;

	/**
	 * Create a new DefaultSourceEdgeDrawing object, given a specific EdgeOntology object.
	 * 
	 * @param teo the TreeEdgeOntology that can semantically interpret the interaction types
	 */
	public DefaultSourceEdgeDrawing(TreeEdgeOntology teo)
	{
		super(teo);
		parentCatToColor = new HashMap<String, Color>();
		parentCatToArrowHead = new HashMap<String, ArrowHead>();
	}

	@Override
	protected Color getEdgeColor(String edgeType)
	{
		if (teo.isDefinedSourceType(edgeType))
		{
			String childCat = teo.getSourceCategory(edgeType);
			Color foundColor = parentCatToColor.get(childCat);
			
			while (foundColor == null && childCat != null)
			{
				String parentCat = teo.retrieveCatParent(childCat);
				foundColor = parentCatToColor.get(parentCat);
				
				childCat = parentCat;
			}
			if (foundColor != null)
			{
				return foundColor;
			}
			return neutral_paint;
		}
		String errormsg = "The provided source type ('" + edgeType + "') is not known in the edge ontology!";
		errormsg += " Make sure to register all source networks to the Project before running the visualisation.";
		throw new IllegalArgumentException(errormsg);
	}

	@Override
	protected ArrowHead getEdgeArrowHead(String edgeType)
	{
		if (teo.isDefinedSourceType(edgeType))
		{
			String childCat = teo.getSourceCategory(edgeType);
			ArrowHead foundArrowHead = parentCatToArrowHead.get(childCat);
			while (foundArrowHead == null && childCat != null)
			{
				String parentCat = teo.retrieveCatParent(childCat);
				foundArrowHead = parentCatToArrowHead.get(parentCat);
				childCat = parentCat;
			}
			if (foundArrowHead != null)
			{
				return foundArrowHead;
			}
			if (teo.isSymmetricalSourceType(edgeType))
			{
				return symm_ah;
			}
			return neutral_ah;
		}
		String errormsg = "The provided source type ('" + edgeType + "') is not known in the edge ontology!";
		errormsg += " Make sure to register all source networks to the Project before running the visualisation.";
		throw new IllegalArgumentException(errormsg);
	}

	/**
	 * Assign a specific paint object to a source category (and its children).
	 * 
	 * @param parentCat a category (also representing its children)
	 * @param p the Color object specifying its visual properties
	 * @throws IllegalArgumentException when the either of the arguments are null, when the type was previously assigned to a paint object, or when the type is not defined in this ontology
	 */
	public void addColor(String parentCat, Color p) throws IllegalArgumentException
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
		if (!teo.isDefinedSourceCat(parentCat))
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
	 * @throws IllegalArgumentException when the either of the arguments are null, when the type was previously assigned to a paint object, or when the type is not defined in this ontology
	 */
	public void addArrowHead(String parentCat, ArrowHead p) throws IllegalArgumentException
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
		if (!teo.isDefinedSourceCat(parentCat))
		{
			String errormsg = "The provided parent category ('" + parentCat + "') is not defined in the edge ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		parentCatToArrowHead.put(parentCat, p);
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
