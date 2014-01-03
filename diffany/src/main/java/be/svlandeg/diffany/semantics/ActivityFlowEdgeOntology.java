package be.svlandeg.diffany.semantics;

import java.awt.Color;
import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.concepts.EdgeDefinition;
import be.svlandeg.diffany.concepts.VisualEdgeStyle.ArrowHead;

/**
 * This edge ontology deals with general flow activity as up/down regulation and
 * their corresponding weights.
 * 
 * @author Sofie Van Landeghem
 */
public class ActivityFlowEdgeOntology extends EdgeOntology
{

	protected String neg_diff_cat;
	protected String pos_diff_cat;

	protected static Color neg_diff_paint = Color.RED;
	protected static Color pos_diff_paint = Color.GREEN;
	protected static Color neutral_diff_paint = Color.GRAY;
	
	protected static ArrowHead neg_diff_ah = ArrowHead.T;
	protected static ArrowHead pos_diff_ah = ArrowHead.ARROW;
	protected static ArrowHead neutral_diff_ah = ArrowHead.ARROW;

	public Set<String> source_pos_cats;
	public Set<String> source_neg_cats;
	public Set<String> source_neutral_cats;

	protected static Color neg_source_paint = Color.RED;
	protected static Color pos_source_paint = Color.GREEN;
	protected static Color neutral_source_paint = Color.LIGHT_GRAY;
	
	protected static ArrowHead neg_source_ah = ArrowHead.T;
	protected static ArrowHead pos_source_ah = ArrowHead.ARROW;
	protected static ArrowHead neutral_source_ah = ArrowHead.ARROW;

	/**
	 * Create a new ontology, defining pos/neg categories. and inserting
	 * After the constructor is called, default edge-category mappings should be
	 * inserted using addCategoryMapping!
	 * 
	 * @param pos_diff_cat the string representing an increase in a differential
	 * network
	 * @param neg_diff_cat the string representing an decrease in a differential
	 * network
	 * @param source_pos_cats positive categories in the input networks (e.g.
	 * upregulation)
	 * @param source_neg_cats negative categories in the input networks (e.g.
	 * inhibition)
	 * @param source_neutral_cats neutral categories in the input networks (e.g.
	 * regulation)
	 */
	public ActivityFlowEdgeOntology(String pos_diff_cat, String neg_diff_cat, Set<String> source_pos_cats, Set<String> source_neg_cats, Set<String> source_neutral_cats)
	{
		super();

		this.pos_diff_cat = pos_diff_cat;
		this.neg_diff_cat = neg_diff_cat;

		addDiffCategory(pos_diff_cat);
		addDiffCategory(neg_diff_cat);

		this.source_pos_cats = source_pos_cats;
		this.source_neg_cats = source_neg_cats;
		this.source_neutral_cats = source_neutral_cats;

		addSourceCategories(source_pos_cats);
		addSourceCategories(source_neg_cats);
		addSourceCategories(source_neutral_cats);
	}

	@Override
	protected Color getDifferentialEdgeColor(String category)
	{
		if (category == null)
		{
			return neutral_diff_paint;
		}
		if (category.equals(pos_diff_cat))
		{
			return pos_diff_paint;
		}
		if (category.equals(neg_diff_cat))
		{
			return neg_diff_paint;
		}
		return neutral_diff_paint;
	}
	
	@Override
	protected ArrowHead getDifferentialEdgeArrowHead(String category)
	{
		if (category == null)
		{
			return neutral_diff_ah;
		}
		if (category.equals(pos_diff_cat))
		{
			return pos_diff_ah;
		}
		if (category.equals(neg_diff_cat))
		{
			return neg_diff_ah;
		}
		return neutral_diff_ah;
	}

	@Override
	protected Color getSourceEdgeColor(String edgeType)
	{
		String cat = getSourceCategory(edgeType);
		if (cat != null && source_pos_cats.contains(cat))
		{
			return pos_source_paint;
		}
		if (cat != null && source_neg_cats.contains(cat))
		{
			return neg_source_paint;
		}
		return neutral_source_paint;
	}

	@Override
	protected ArrowHead getSourceEdgeArrowHead(String edgeType)
	{
		String cat = getSourceCategory(edgeType);
		if (cat != null && source_pos_cats.contains(cat))
		{
			return pos_source_ah;
		}
		if (cat != null && source_neg_cats.contains(cat))
		{
			return neg_source_ah;
		}
		return neutral_source_ah;
	}

	@Override
	public EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, Set<EdgeDefinition> conEdges, double cutoff) throws IllegalArgumentException
	{
		// Note: this implementation could also be done by refining the algorithm of ProcessOntology.
		// The current implementation does the same in principle, but basically
		// interprets positive/negative edges straight away instead of at the end.

		EdgeDefinition diff_edge = new EdgeDefinition();
		Set<EdgeDefinition> conEdges2 = new HashSet<EdgeDefinition>();

		// ////////// DEAL WITH SYMMETRY AND NEGATION ////////////////
		boolean conSymm = true;
		for (EdgeDefinition e : conEdges)
		{
			// negated edges are set to void
			if (e.isNegated())
			{
				conEdges2.add(EdgeDefinition.getVoidEdge());
			} else
			{
				conEdges2.add(e);
			}
			if (!e.isSymmetrical())
			{
				// the set of condition-specific edges is only symmetrical when all edges are
				conSymm = false;
			}
		}

		if (conEdges.isEmpty())
		{
			conEdges2.add(EdgeDefinition.getVoidEdge());
		}

		boolean refNeg = refEdge.isNegated();
		if (refNeg)
			refEdge = EdgeDefinition.getVoidEdge();

		diff_edge.makeNegated(false); // a differential edge is never negated

		boolean refSymm = refEdge.isSymmetrical();
		boolean diffSymm = refSymm && conSymm;
		diff_edge.makeSymmetrical(diffSymm);

		// ////////// DETERMINE TYPE AND WEIGHT ////////////////

		String refCat = getSourceCategory(refEdge.getType());
		boolean refVoid = refCat.equals(VOID_TYPE);
		// the differential edge can not be calculated if the reference edge is
		// a neutral edge (unspecified)
		if (source_neutral_cats.contains(refCat))
		{
			return EdgeDefinition.getVoidEdge();
		}

		Boolean refDirection = null;
		if (source_pos_cats.contains(refCat))
			refDirection = true;
		else if (source_neg_cats.contains(refCat))
			refDirection = false;

		int countUp = 0;
		int countDown = 0;
		int countNeutral = 0;
		double minDiffWeight = Double.MAX_VALUE;

		for (EdgeDefinition conEdge : conEdges2)
		{
			String conCat = getSourceCategory(conEdge.getType());
			double diffWeight = conEdge.getWeight() - refEdge.getWeight();

			if (source_pos_cats.contains(conCat))
			{
				if (refVoid)
				{
					diffWeight = conEdge.getWeight();
					countUp++;
				} else
				{
					if (refDirection && diffWeight > 0)
					{
						countUp++;
					}
					if (refDirection && diffWeight < 0)
					{
						countDown++;
					}
					if (!refDirection)
					{
						countUp++;
						diffWeight = conEdge.getWeight() + refEdge.getWeight();
					}
				}
			} else if (source_neg_cats.contains(conCat))
			{
				if (refVoid)
				{
					countDown++;
					diffWeight = conEdge.getWeight();
				} else
				{
					if (!refDirection && diffWeight > 0)
					{
						countDown++;
					}
					if (!refDirection && diffWeight < 0)
					{
						countUp++;
					}
					if (refDirection)
					{
						countDown++;
						diffWeight = conEdge.getWeight() + refEdge.getWeight();
					}
				}
			} else if (conCat.equals(VOID_TYPE))
			{
				diffWeight = refEdge.getWeight();
				if (refVoid)
				{
					countNeutral++;
				} else if (refDirection)
				{
					countDown++;
				} else
				{
					countUp++;
				}
			} else if (source_neutral_cats.contains(conCat))
			{
				countNeutral++;
			}
			minDiffWeight = Math.min(minDiffWeight, diffWeight);
		}

		// some are positive, some are negative -> no general differential edge
		if (countUp > 0 && countDown > 0)
		{
			return EdgeDefinition.getVoidEdge();
		}

		Boolean conDirection = null;
		if (countDown == conEdges2.size())
		{
			// all differential edges are negative (or void)
			conDirection = false;
		}
		if (countUp == conEdges2.size())
		{
			// all differential edges are positive (or void)
			conDirection = true;
		}

		// the differential edge can not be calculated if there are neutral/unspecified edges
		if (conDirection == null || countNeutral > 0)
		{
			return EdgeDefinition.getVoidEdge();
		}

		if (minDiffWeight <= cutoff)
		{
			return EdgeDefinition.getVoidEdge();
		}

		if (conDirection)
			diff_edge.setType(pos_diff_cat);
		if (!conDirection)
			diff_edge.setType(neg_diff_cat);

		diff_edge.setWeight(minDiffWeight);

		return diff_edge;
	}

}
