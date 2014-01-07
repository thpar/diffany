package be.svlandeg.diffany.concepts;

import java.util.*;


/**
 * This class is a simplification of EdgeSet, containing only a single edge between two nodes in any of the networks.
 * It knows which edges are from the reference network and from the condition-dependent networks,
 * and it knows the number of condition-specific networks.
 * 
 * @author Sofie Van Landeghem
 */
public class SingleEdgeSet
{
	
	private EdgeDefinition referenceEdge;
	private List<EdgeDefinition> conditionEdges;
	private int conditionCount;
	
	/**
	 * Constructor: generates a new SingleEdgeSet
	 * @param referenceEdge the reference edges
	 * @param conditionEdges the condition-specific edges, 1 (may be void) for each condition
	 */
	public SingleEdgeSet(Edge referenceEdge, List<Edge> conditionEdges)
	{
		this(conditionEdges.size());
		this.referenceEdge = referenceEdge;
		int i = 0; 
		for (Edge e : conditionEdges)
		{
			putConditionEdge(e, i);
			i++;
		}
	}
	
	/**
	 * Constructor: generates an empty SingleEdgeSet
	 * @param conditionCount the number of condition-specific networks
	 */
	public SingleEdgeSet(int conditionCount)
	{
		referenceEdge = null;
		conditionEdges = new ArrayList<EdgeDefinition>();
		for (int i = 0; i < conditionCount; i++)
		{
			conditionEdges.add(i, null);
		}
		this.conditionCount = conditionCount;
	}
	
	/**
	 * Return the reference edge in this set
	 * @return the reference edge in this set
	 */
	public EdgeDefinition getReferenceEdge()
	{
		return referenceEdge;
	}

	/**
	 * Return all condition-specific edges in this set
	 * @return all condition-specific edges in this set
	 */
	public List<EdgeDefinition> getConditionEdges()
	{
		return conditionEdges;
	}
	
	/**
	 * Return the condition-specific edge for one specific condition-dependent network in this set
	 * @return the correct condition-specific edge
	 */
	public EdgeDefinition getConditionEdge(int conditionNR)
	{
		if (conditionNR < 0 || conditionNR >= getConditionCount())
		{
			throw new IllegalArgumentException("The condition number should be larger than 0 and smaller than the number of conditions " + conditionNR + "!");
		}
		return conditionEdges.get(conditionNR);
	}
	
	/**
	 * Return all edges in this set: both reference as well as condition-specific edges
	 * @return all edges in this set
	 */
	public Set<EdgeDefinition> getAllEdges()
	{
		Set<EdgeDefinition> allEdges = new HashSet<EdgeDefinition>();
		for (EdgeDefinition s : conditionEdges)
		{
			allEdges.add(s);
		}
		allEdges.add(referenceEdge);
		return allEdges;
	}

	/**
	 * Get the number of condition-specific networks 
	 * @return the number of condition-specific networks 
	 */
	public int getConditionCount()
	{
		return conditionCount;
	}
	
	/**
	 * Define the reference edge
	 * @param the new reference edge
	 */
	public void putReferenceEdge(EdgeDefinition e)
	{
		referenceEdge = e;
	}
	
	/**
	 * Add a condition-specific edge
	 * @param e a new condition-specific edge
	 * @param conditionNR the number of the condition-specific network
	 * @throws 
	 */
	public void putConditionEdge(EdgeDefinition e, int conditionNR)
	{
		if (conditionNR < 0 || conditionNR >= getConditionCount())
		{
			throw new IllegalArgumentException("The condition number should be larger than 0 and smaller than the number of conditions !");
		}
		conditionEdges.set(conditionNR,e);
	}

	@Override
	public String toString()
	{
		String result = "Single Edge Set with reference edge: [" + referenceEdge + "] ";
		result += " and " + conditionCount + " conditions: ";
		for (EdgeDefinition e : conditionEdges)
		{
			result += "[" + e + "] ";
		}
		return result;
	}
	

}
