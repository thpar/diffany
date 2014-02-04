package be.svlandeg.diffany.concepts;

import java.util.*;

/**
 * This class represents an edge set between two previously defined nodes: all edges between those nodes. It knows which edges are from the reference network and from the condition-dependent networks, and it knows the number of condition-specific networks.
 * 
 * @author Sofie Van Landeghem
 */
public class EdgeSet
{

	private Set<EdgeDefinition> referenceEdges;
	private List<Set<EdgeDefinition>> conditionEdges;
	private int conditionCount;

	/**
	 * Constructor: generates a new EdgeSet
	 * 
	 * @param referenceEdges the reference edges
	 * @param conditionEdges the condition-specific edges, 1 set (may be empty) for each condition
	 */
	public EdgeSet(Set<Edge> referenceEdges, List<Set<Edge>> conditionEdges)
	{
		this(conditionEdges.size());
		for (Edge e : referenceEdges)
		{
			addReferenceEdge(e);
		}
		for (int i = 0; i < conditionEdges.size(); i++)
		{
			for (Edge e : conditionEdges.get(i))
			{
				addConditionEdge(e, i);
			}
		}
	}

	/**
	 * Constructor: generates an empty EdgeSet
	 * 
	 * @param conditionCount the number of condition-specific networks
	 */
	public EdgeSet(int conditionCount)
	{
		this.conditionCount = conditionCount;
		referenceEdges = new HashSet<EdgeDefinition>();
		conditionEdges = new ArrayList<Set<EdgeDefinition>>();
		for (int i = 0; i < conditionCount; i++)
		{
			conditionEdges.add(i, new HashSet<EdgeDefinition>());
		}

	}

	/**
	 * Return all reference edges in this set
	 * 
	 * @return all reference edges in this set
	 */
	public Set<EdgeDefinition> getReferenceEdges()
	{
		return referenceEdges;
	}

	/**
	 * Return all condition-specific edges in this set
	 * 
	 * @return all condition-specific edges in this set
	 */
	public List<Set<EdgeDefinition>> getConditionEdges()
	{
		return conditionEdges;
	}

	/**
	 * Return all condition-specific edges for one specific condition-dependent network in this set
	 * 
	 * @return all condition-specific edges
	 */
	public Set<EdgeDefinition> getConditionEdges(int conditionNR)
	{
		if (conditionNR < 0 || conditionNR >= getConditionCount())
		{
			throw new IllegalArgumentException("The condition number should be larger than 0 and smaller than the number of conditions " + conditionNR + "!");
		}
		return conditionEdges.get(conditionNR);
	}

	/**
	 * Return all edges in this set: both reference as well as condition-specific edges
	 * 
	 * @return all edges in this set
	 */
	public Set<EdgeDefinition> getAllEdges()
	{
		Set<EdgeDefinition> allEdges = new HashSet<EdgeDefinition>();
		for (Set<EdgeDefinition> s : conditionEdges)
		{
			allEdges.addAll(s);
		}
		allEdges.addAll(referenceEdges);
		return allEdges;
	}

	/**
	 * Get the number of condition-specific networks
	 * 
	 * @return the number of condition-specific networks
	 */
	public int getConditionCount()
	{
		return conditionCount;
	}

	/**
	 * Add a reference edge
	 * 
	 * @param e a new reference edge
	 */
	public void addReferenceEdge(EdgeDefinition e)
	{
		referenceEdges.add(e);
	}

	/**
	 * Add a condition-specific edge
	 * 
	 * @param e a new condition-specific edge
	 * @param conditionNR the number of the condition-specific network
	 * @throws IllegalArgumentException when the number of the condition-specific network is invalid
	 */
	public void addConditionEdge(EdgeDefinition e, int conditionNR)
	{
		if (conditionNR < 0 || conditionNR >= getConditionCount())
		{
			throw new IllegalArgumentException("The condition number should be larger than 0 and smaller than the number of conditions !");
		}
		conditionEdges.get(conditionNR).add(e);
	}

	@Override
	public String toString()
	{
		String result = "Edge Set with reference edges: [";
		for (EdgeDefinition e : referenceEdges)
		{
			result += e + " ";
		}
		result += "] ";
		result += " and " + conditionCount + " conditions: ";
		for (int i = 0; i < conditionCount; i++)
		{
			result += " [";
			for (EdgeDefinition e : conditionEdges.get(i))
			{
				result += e + " ";
			}
			result += "] ";
		}
		return result;
	}

}
