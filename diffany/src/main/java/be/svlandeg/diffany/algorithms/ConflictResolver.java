package be.svlandeg.diffany.algorithms;

import java.util.*;

import be.svlandeg.diffany.concepts.EdgeDefinition;
import be.svlandeg.diffany.concepts.EdgeSet;
import be.svlandeg.diffany.concepts.SingleEdgeSet;
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
	 * TODO
	 * @param eo
	 * @param oldEdgeSet
	 * @return
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

}
