package be.svlandeg.diffany.algorithms;

import java.util.*;

import be.svlandeg.diffany.concepts.*;
import be.svlandeg.diffany.semantics.EdgeOntology;

/**
 * This class provides methods to resolve conflicts in the input networks, 
 * e.g. positive and negative regulation between two nodes,
 * or PTM and phosphorylation at the same time.
 * 
 * @author Sofie Van Landeghem
 */
public class ConflictResolver
{
	
	/**
	 * Group all input edges into subclasses per root category of the EdgeOntology
	 * @param eo the edge ontology
	 * @param oldEdgeSet the original sets of input edges
	 * @return all input edges grouped by edge root category
	 */
	public Map<String, SingleEdgeSet> resolveEdgesPerRoot(EdgeOntology eo, EdgeSet oldEdgeSet)
	{
		Map<String, SingleEdgeSet> mappedEdges = new HashMap<String, SingleEdgeSet>();
		Set<String> roots = eo.retrieveAllSourceRootCats();

		Set<EdgeDefinition> refEs = oldEdgeSet.getReferenceEdges();
		List<Set<EdgeDefinition>> conditionEdges = oldEdgeSet.getConditionEdges();
		
		for (String root : roots)
		{
			SingleEdgeSet es = new SingleEdgeSet(oldEdgeSet.getConditionCount());
			
			// Determine which reference edges belong to this root
			Set<EdgeDefinition> rootRefEs = new HashSet<EdgeDefinition>();
			for (EdgeDefinition refE : refEs)
			{
				String edgeType = refE.getType();
				String edgeClass = eo.getSourceCategory(edgeType);
				boolean belongsToRoot = eo.isSourceChildOf(edgeClass, root) >= 0;
				if (belongsToRoot)
				{
					rootRefEs.add(refE);
				}
			}
			if (rootRefEs.size() > 1)	// TODO
			{
				throw new UnsupportedOperationException("This algorithm currently only supports 1 edge between two nodes in the original networks");
			}
			if (rootRefEs.size() < 1)
			{
				es.putReferenceEdge(EdgeDefinition.getVoidEdge());
			}
			else
			{
				es.putReferenceEdge(rootRefEs.iterator().next());
			}
			
			// Determine which condition-specific edges belong to this root
			for (int i = 0; i < oldEdgeSet.getConditionCount(); i++)
			{
				Set<EdgeDefinition> rootCondIEs = new HashSet<EdgeDefinition>();
				for (EdgeDefinition conE : conditionEdges.get(i))
				{
					String edgeType = conE.getType();
					String edgeClass = eo.getSourceCategory(edgeType);
					boolean belongsToRoot = eo.isSourceChildOf(edgeClass, root) >= 0;
					if (belongsToRoot)
					{
						rootCondIEs.add(conE);
					}
				}
				if (rootCondIEs.size() > 1)	// TODO
				{
					throw new UnsupportedOperationException("This algorithm currently only supports 1 edge between two nodes in the original networks");
				}
				if (rootCondIEs.size() < 1)
				{
					es.putConditionEdge(EdgeDefinition.getVoidEdge(), i);
				}
				else
				{
					es.putConditionEdge(rootCondIEs.iterator().next(), i);
				}
				
			}
			mappedEdges.put(root, es);
		}
		
		return mappedEdges;
	}
	
	// TODO
	public SingleEdgeSet unifyDirection(SingleEdgeSet oldSet, EdgeOntology eo)
	{
		int conditions = oldSet.getConditionCount();
		
		EdgeDefinition old_referenceEdge = oldSet.getReferenceEdge();
		List<EdgeDefinition> old_conditionEdges = oldSet.getConditionEdges();
		
		int symmetricalCount = 0;
		int directedCount = 0;
		
		if (old_referenceEdge.isSymmetrical())
		{
			symmetricalCount++;
		}
		else
		{
			directedCount++;
		}
		
		for (EdgeDefinition c : old_conditionEdges)
		{
			if (c.isSymmetrical())
			{
				symmetricalCount++;
			}
			else
			{
				directedCount++;
			}
		}
		
		if (symmetricalCount == 0 || directedCount == 0)
		{
			return oldSet;	// nothing needs to be changed
		}
		
		// let's make all of them directed!
		SingleEdgeSet newSet = new SingleEdgeSet(conditions);
		
		String symmType = old_referenceEdge.getType();

		EdgeDefinition fwdDirection = new EdgeDefinition(old_referenceEdge);
		fwdDirection.makeSymmetrical(false);
		
		newSet.putReferenceEdge(fwdDirection);
		
		return newSet;
		
	}
	

}
