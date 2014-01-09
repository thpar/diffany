package be.svlandeg.diffany.algorithms;

import java.util.*;

import be.svlandeg.diffany.concepts.*;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * This class provides generic methods useful for network cleaning before or
 * after applying differential algorithms.
 * 
 * @author Sofie Van Landeghem
 */
public class NetworkCleaning
{
	
	private Logger log;
	
	/**
	 * Create a new cleaning object, which can log important messages.
	 * @param log the logger object
	 */
	public NetworkCleaning(Logger log)
	{
		this.log = log;
	}
	
	/**
	 * Clean an output network: Remove redundant/duplicate edges in the network.
	 * 
	 * Calls the subroutine removeRedundantEdges.
	 * 
	 * @param net the network that needs cleaning
	 * @param nm the node mapper for the network
	 **/
	public void fullOutputCleaning(Network net, NodeMapper nm)
	{
		removeRedundantEdges(net, nm);
	}

	/**
	 * Remove edges in the network that are symmetrical and are represented
	 * twice (source-target and target-source). One of the two is removed only
	 * then when the type, weight and negation are all equal.
	 * 
	 * @param net the network that needs cleaning
	 * @param nm the node mapper for the network
	 **/
	private void removeRedundantEdges(Network net, NodeMapper nm)
	{
		Set<Edge> removed_edges = new HashSet<Edge>();
		
		// remove duplicate symmetrical edges between source-target and
		// target-source
		for (Node n1 : net.getNodes())
		{
			String name1 = n1.getName();
			for (Node n2 : net.getNodes())
			{
				String name2 = n2.getName();
				if (name1.compareTo(name2) < 0)
				{
					Set<Edge> all_edges = net.getAllEdges(n1, n2, nm);
					for (Edge et : all_edges)
					{
						for (Edge eb : all_edges)
						{
							if (!et.equals(eb) && !removed_edges.contains(et))
							{
								if ((et.isSymmetrical() && eb.isSymmetrical()) && (et.getType().equals(eb.getType())))
								{
									if ((et.getWeight() == eb.getWeight()) && (et.isNegated() == eb.isNegated()))
									{
										net.removeEdge(eb);
										removed_edges.add(eb);
									}
								}
							}
						}
					}
				}
			}
		}
		log.log("Output cleaning: removed redundant edges");
	}
	
	/**
	 * Clean an input network:
	 * Group all input edges into subclasses per root category of the EdgeOntology, unify the directionality 
	 * (either all symmetric or all directed), and resolve conflicts within a root category.
	 * 
	 * Calls the subroutines resolveEdgesPerRoot, unifyDirection and summarizeToOne.
	 * 
	 * @param eo the edge ontology
	 * @param oldEdgeSet all edges between two nodes, including both directed and symmetrical edges
	 * @param backEdgeSet all reverse edges between the same two nodes, excluding symmetrical edges
	 * @return all edges grouped by semantic root category, with unified directionality, and only 1 edge per network and root cat.
	 */
	public Map<String, SingleEdgeSet> fullInputCleaning(EdgeOntology eo, EdgeSet oldEdgeSet, EdgeSet backEdgeSet)
	{
		Map<String, EdgeSet> mappedNormalEdges = resolveEdgesPerRoot(eo, oldEdgeSet);
		Map<String, EdgeSet> mappedReverseEdges = resolveEdgesPerRoot(eo, backEdgeSet);
		
		Map<String, SingleEdgeSet> mappedSingleEdges = new HashMap<String, SingleEdgeSet>();
		
		for (String key : mappedNormalEdges.keySet())
		{
			EdgeSet normalEdges = mappedNormalEdges.get(key);
			EdgeSet reverseEdges = mappedReverseEdges.get(key);
			
			EdgeSet unifiedEdgeSet = unifyDirection(normalEdges, reverseEdges, eo);
			
			SingleEdgeSet singleEdgeSet = summarizeToOne(unifiedEdgeSet,eo);
			
			mappedSingleEdges.put(key, singleEdgeSet);
		}
		return mappedSingleEdges;
	}
	

	/**
	 * Group all input edges into subclasses per root category of the EdgeOntology.
	 *
	 * @param eo the edge ontology
	 * @param oldEdgeSet the original sets of input edges
	 * @return all input edges grouped by edge root category
	 */
	protected Map<String, EdgeSet> resolveEdgesPerRoot(EdgeOntology eo, EdgeSet oldEdgeSet)
	{
		Map<String, EdgeSet> mappedEdges = new HashMap<String, EdgeSet>();
		Set<String> roots = eo.retrieveAllSourceRootCats();

		Set<EdgeDefinition> refEs = oldEdgeSet.getReferenceEdges();
		List<Set<EdgeDefinition>> conditionEdges = oldEdgeSet.getConditionEdges();
		
		for (String root : roots)
		{
			EdgeSet es = new EdgeSet(oldEdgeSet.getConditionCount());
			
			// Determine which reference edges belong to this root
			for (EdgeDefinition refE : refEs)
			{
				String edgeType = refE.getType();
				String edgeClass = eo.getSourceCategory(edgeType);
				boolean belongsToRoot = eo.isSourceChildOf(edgeClass, root) >= 0;
				if (belongsToRoot)
				{
					es.addReferenceEdge(refE);
				}
			}
			if (es.getReferenceEdges().isEmpty())
			{
				es.addReferenceEdge(EdgeDefinition.getVoidEdge());
			}
			
			// Determine which condition-specific edges belong to this root
			for (int i = 0; i < oldEdgeSet.getConditionCount(); i++)
			{
				for (EdgeDefinition conE : conditionEdges.get(i))
				{
					String edgeType = conE.getType();
					String edgeClass = eo.getSourceCategory(edgeType);
					boolean belongsToRoot = eo.isSourceChildOf(edgeClass, root) >= 0;
					if (belongsToRoot)
					{
						es.addConditionEdge(conE, i);
					}
				}
				if (es.getConditionEdges(i).isEmpty())
				{
					es.addConditionEdge(EdgeDefinition.getVoidEdge(), i);
				}
			}
			mappedEdges.put(root, es);
		}
		
		return mappedEdges;
	}
	
	
	/**
	 * Create a new EdgeSet with all interactions either directed, or symmetrical.
	 * This should only be done for EdgeSets within the same ontology subtree (i.e. call resolveEdgesPerRoot first).
	 * 
	 * @param oldSet the old edgeset which might have a mixture of directed and symmetrical edges
	 * @param oldSet the old edgeset with opposite directed edges (if any)
	 * @param eo the edge ontology
	 * @return the new edge set with only directed, or only symmetrical edges
	 */
	protected EdgeSet unifyDirection(EdgeSet oldSet, EdgeSet backEdgeSet, EdgeOntology eo)
	{
		boolean hasSymmetrical = containsSymmetricalEdges(oldSet);
		boolean hasDirected = containsDirectedEdges(oldSet);
		boolean hasOppositeDirected = containsDirectedEdges(backEdgeSet);
		
		// all edges are directed
		if (! hasSymmetrical)
		{
			return oldSet;	// nothing needs to be changed
		}
		
		// all edges are symmetrical, also the reverse ones
		if (! hasDirected && ! hasOppositeDirected)
		{
			return oldSet;	// nothing needs to be changed
		}
		
		// In all other cases: let's make all of them directed!
		int conditions = oldSet.getConditionCount();
		
		EdgeSet newSet = new EdgeSet(conditions);
		
		for (EdgeDefinition referenceEdge : oldSet.getReferenceEdges())
		{
			EdgeDefinition refFwdDirection = new EdgeDefinition(referenceEdge);
			refFwdDirection.makeSymmetrical(false);
			newSet.addReferenceEdge(refFwdDirection);
		}
		
		for (int i = 0; i < conditions; i++)
		{
			for (EdgeDefinition old_conIEdge : oldSet.getConditionEdges(i))
			{
				EdgeDefinition conIFwdDirection = new EdgeDefinition(old_conIEdge);
				conIFwdDirection.makeSymmetrical(false);
				newSet.addConditionEdge(conIFwdDirection, i);
			}
		}
		
		return newSet;
	}
	
	/**
	 * Summarize all edges per network to one edge only.
	 * It is assumed that unifyDirection was previously called to obtain either all directed or all symmetrical edges,
	 * and that resolveEdgesPerRoot was used to provide an EdgeSet which only contains edges for one root category.
	 * 
	 * @param oldSet the old edge set
	 * @param eo the edge ontology
	 * @return the new edge set, holding at most one edge per input network
	 */
	protected SingleEdgeSet summarizeToOne(EdgeSet oldSet, EdgeOntology eo)
	{
		SingleEdgeSet newSet = new SingleEdgeSet(oldSet.getConditionCount());
		
		Set<EdgeDefinition> rootRefEs = oldSet.getReferenceEdges();
		if (rootRefEs.size() > 1)	// TODO
		{
			newSet.putReferenceEdge(resolveToOne(rootRefEs,  eo));
		}
		else
		{
			newSet.putReferenceEdge(rootRefEs.iterator().next());
		}
		
		for (int i = 0; i < oldSet.getConditionCount(); i++)
		{
			Set<EdgeDefinition> rootCondIEs = oldSet.getConditionEdges(i);
			if (rootCondIEs.size() > 1)	// TODO
			{
				newSet.putConditionEdge(resolveToOne(rootCondIEs, eo), i);
			}
			else
			{
				newSet.putConditionEdge(rootCondIEs.iterator().next(), i);
			}
		}
		
		return newSet;
	}
	
	/**
	 * Resolve a set of edges to one. This is currently implemented by taking the edge with the highest weight. 
	 * 
	 * @param edges the original set of input edges
	 * @param eo the edge ontology
	 * @return one edge, produced after resolving conflicts, or throws a RunTimeException if no best edge could be found. 
	 */
	protected EdgeDefinition resolveToOne(Set<EdgeDefinition> edges,  EdgeOntology eo)
	{
		// TODO: should we also take into account whether or one of the edges is more specific?
		double maxWeight = 0.0;
		
		for (EdgeDefinition e : edges)
		{
			maxWeight = Math.max(maxWeight, e.getWeight());
		}
		
		for (EdgeDefinition e : edges)
		{
			if (maxWeight == e.getWeight())
			{
				return e;
			}
		}
		String errormsg = "Could not resolve the set of edges to one!";
		throw new RuntimeException(errormsg);
	}
	
	/**
	 * Determine whether or not the set of edges contains at least one directed edge
	 * @param edgeSet the set of edges
	 * @return whether or not the set contains at least one directed edge
	 */
	private boolean containsDirectedEdges(EdgeSet edgeSet)
	{
		int directedCount = 0;
		
		for (EdgeDefinition referenceEdge : edgeSet.getReferenceEdges())
		{
			if (! referenceEdge.isSymmetrical())
			{
				directedCount++;
			}
		}
		
		List<Set<EdgeDefinition>> conditionEdges = edgeSet.getConditionEdges();
		for (int i = 0; i < conditionEdges.size(); i++)
		{
			for (EdgeDefinition c : conditionEdges.get(i))
			{
				if (! c.isSymmetrical())
				{
					directedCount++;
				}
			}
		}
		return (directedCount > 0);
	}
	
	/**
	 * Determine whether or not the set of edges contains at least one symmetrical edge
	 * @param edgeSet the set of edges
	 * @return whether or not the set contains at least one symmetrical edge
	 */
	private boolean containsSymmetricalEdges(EdgeSet edgeSet)
	{
		int symmetricalCount = 0;
		
		for (EdgeDefinition referenceEdge : edgeSet.getReferenceEdges())
		{
			if (referenceEdge.isSymmetrical())
			{
				symmetricalCount++;
			}
		}
		
		List<Set<EdgeDefinition>> conditionEdges = edgeSet.getConditionEdges();
		for (int i = 0; i < conditionEdges.size(); i++)
		{
			for (EdgeDefinition c : conditionEdges.get(i))
			{
				if (c.isSymmetrical())
				{
					symmetricalCount++;
				}
			}
		}
		return (symmetricalCount > 0);
	}
	

}
