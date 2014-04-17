package be.svlandeg.diffany.core.networks.merged;

import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.core.networks.Edge;

/**
 * This class allows quick conversion between merged input networks with condition edges, and normal input networks 
 * with normal edges.
 * 
 * @author Sofie Van Landeghem
 */
public class MergedConvertor
{

	/**
	 * Convert a set of condition edges to a set of normal edges.
	 * @param cEdges the condition-specific edges
	 * @return the set of normal edges
	 */
	public static Set<Edge> convertToNormalEdges(Set<ConditionEdge> cEdges)
	{
		Set<Edge> edges = new HashSet<Edge>();
		for (ConditionEdge e : cEdges)
		{
			edges.add(e);
		}
		return edges;
	}
	
	/**
	 * Convert a set of normal edges to a set of condition edges.
	 * @param edges a set of normal edges which are actually condition-specific (and can thus be cast!)
	 * @return the set of condition-specific edges
	 * @throws ClassCastException when the set of input edges are not of the type ConditionEdge
	 */
	public static Set<ConditionEdge> convertToConditionEdges(Set<Edge> edges)
	{
		Set<ConditionEdge> cEdges = new HashSet<ConditionEdge>();
		for (Edge e : edges)
		{
			cEdges.add((ConditionEdge) e);
		}
		return cEdges;
	}

}
