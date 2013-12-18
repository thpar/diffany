package be.svlandeg.diffany.algorithms;

import java.util.*;

import be.svlandeg.diffany.concepts.*;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * This class can calculate a differential network between one reference and one condition-specific network. 
 * Currently this algorithm assumes a 1-to-1 mapping of nodes between the two networks!
 * 
 * @author Sofie Van Landeghem
 */
public class CalculateDiffOfTwo
{
	
	protected NetworkCleaning cleaning;
	
	/**
	 * Constructor, which initializes the functionality for cleaning a network.
	 */
	public CalculateDiffOfTwo()
	{
		cleaning = new NetworkCleaning();
	}
	
	
	/**
	 * Calculate the differential network between the reference and condition-specific network. 
	 * The overlapping network should be calculated independently!
	 * This method can only be called from within the package (CalculateDiff) and can thus assume proper input.
	 * 
	 * @param reference the reference network
	 * @param condition a condition-specific network
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param diff_name the name to give to the differential network. 
	 * @return the differential network between the two
	 *         
	 * TODO: expand this algorithm to be able to deal with n-m node mappings (v.2.0) 
	 * TODO: expand this algorithm to be able to deal with more than 1 edge between two nodes 
	 * in the original networks (v.1.0)
	 */
	protected DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, ConditionNetwork condition, EdgeOntology eo, 
			NodeMapper nm, String diff_name, double cutoff)
	{
		Set<ConditionNetwork> conditions = new HashSet<ConditionNetwork>();
		conditions.add(condition);

		Set<Network> allOriginals = new HashSet<Network>();
		allOriginals.addAll(conditions);
		allOriginals.add(reference);

		DifferentialNetwork diff = new DifferentialNetwork(diff_name, reference, conditions);

		Map<Node, Set<Node>> nodeMapping = nm.getAllEquals(reference, condition);
		Set<Node> allNodes = nm.getAllNodes(allOriginals);

		Map<String, Node> allDiffNodes = new HashMap<String, Node>();

		for (Node source1 : allNodes) // source node in reference network
		{
			// get the equivalent source node in the condition network
			Node source2;
			if (nodeMapping.containsKey(source1))
			{
				Set<Node> sources2 = nodeMapping.get(source1);
				source2 = getSingleNode(sources2);
			}
			else
			// source1 is not actually a part of the reference network
			{
				source2 = source1;
			}

			for (Node target1 : allNodes) // target node in reference network
			{
				// get the equivalent target node in the condition network
				Node target2;
				if (nodeMapping.containsKey(target1))
				{
					Set<Node> targets2 = nodeMapping.get(target1);
					target2 = getSingleNode(targets2);
				}
				else
				// target1 is not actually a part of the reference network
				{
					target2 = target1;
				}

				// get the reference edge
				Set<Edge> referenceEdges = reference.getAllEdges(source1, target1);
				EdgeDefinition edgedef1 = getSingleEdge(referenceEdges);

				// get the condition-specific edge
				Set<Edge> conditionEdges = new HashSet<Edge>();
				if (source2 != null && target2 != null)
				{
					conditionEdges = condition.getAllEdges(source2, target2);
				}
				EdgeDefinition edgedef2 = getSingleEdge(conditionEdges);
				Set<EdgeDefinition> edgedefs2 = new HashSet<EdgeDefinition>();
				edgedefs2.add(edgedef2);

				EdgeDefinition diff_edge_def = eo.getDifferentialEdge(edgedef1, edgedefs2, cutoff);

				Set<Node> allSources = new HashSet<Node>();
				allSources.add(source1);
				allSources.add(source2);
				String sourceconsensus = nm.getConsensusName(allSources);
				
				Set<Node> allTargets = new HashSet<Node>();
				allTargets.add(target1);
				allTargets.add(target2);
				String targetconsensus = nm.getConsensusName(allTargets);
				
				// non-void differential edge
				if (diff_edge_def.getType() != EdgeOntology.VOID_TYPE && sourceconsensus != null && targetconsensus != null)
				{
					if (!allDiffNodes.containsKey(sourceconsensus))
					{
						allDiffNodes.put(sourceconsensus, new Node(sourceconsensus));
					}
					Node sourceresult = allDiffNodes.get(sourceconsensus);
					
					if (!allDiffNodes.containsKey(targetconsensus))
					{
						allDiffNodes.put(targetconsensus, new Node(targetconsensus));
					}
					Node targetresult = allDiffNodes.get(targetconsensus);
					
					Edge edgediff = new Edge(sourceresult, targetresult, diff_edge_def);
					diff.addEdge(edgediff);
				}
			}
		}
		cleaning.removeRedundantEdges(diff);
		cleaning.directSymmetricalWhenOverlapping(diff);
		
		return diff;
	}
	
	/**
	 * Calculate the overlapping network between two networks.
	 * This method can only be called from within the package (CalculateDiff) and can thus assume proper input.
	 * 
	 * @param n1 the first network
	 * @param n2 the second network
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param overlap_name the name to give to the overlapping network
	 * @param minOperator whether or not to take the minimum of the edge weights - if false, the maximum is taken
	 * 
	 * @return the overlapping network between the two
	 *      
	 * TODO: expand this algorithm to be able to deal with n-m node mappings (v.2.0) 
	 * TODO: expand this algorithm to be able to deal with more than 1 edge between two nodes 
	 * in the original networks (v.1.0)
	 */
	protected OverlappingNetwork calculateOverlappingNetwork(Network n1, Network n2, EdgeOntology eo, 
			NodeMapper nm, String overlap_name, double cutoff, boolean minOperator)
	{
		Set<Network> allOriginals = new HashSet<Network>();
		allOriginals.add(n1);
		allOriginals.add(n2);

		OverlappingNetwork overlap = new OverlappingNetwork(overlap_name, allOriginals);

		Map<Node, Set<Node>> nodeMapping = nm.getAllEquals(n1, n2);
		Set<Node> allNodes = nm.getAllNodes(allOriginals);

		Map<String, Node> allDiffNodes = new HashMap<String, Node>();

		for (Node source1 : allNodes) // source node in reference network
		{
			// get the equivalent source node in the condition network
			Node source2;
			if (nodeMapping.containsKey(source1))
			{
				Set<Node> sources2 = nodeMapping.get(source1);
				source2 = getSingleNode(sources2);
			}
			else
			// source1 is not actually a part of the reference network
			{
				source2 = source1;
			}

			for (Node target1 : allNodes) // target node in reference network
			{
				// get the equivalent target node in the condition network
				Node target2;
				if (nodeMapping.containsKey(target1))
				{
					Set<Node> targets2 = nodeMapping.get(target1);
					target2 = getSingleNode(targets2);
				}
				else
				// target1 is not actually a part of the reference network
				{
					target2 = target1;
				}

				// get the reference edge
				Set<Edge> n1Edges = n1.getAllEdges(source1, target1);
				EdgeDefinition edgedef1 = getSingleEdge(n1Edges);

				// get the condition-specific edge
				Set<Edge> n2Edges = new HashSet<Edge>();
				if (source2 != null && target2 != null)
				{
					n2Edges = n2.getAllEdges(source2, target2);
				}
				EdgeDefinition edgedef2 = getSingleEdge(n2Edges);

				Set<EdgeDefinition> edges = new HashSet<EdgeDefinition>();
				edges.add(edgedef1);
				edges.add(edgedef2);
				EdgeDefinition overlap_edge_def = eo.getOverlapEdge(edges, cutoff, minOperator);

				Set<Node> allSources = new HashSet<Node>();
				allSources.add(source1);
				allSources.add(source2);
				String sourceconsensus = nm.getConsensusName(allSources);
				
				Set<Node> allTargets = new HashSet<Node>();
				allTargets.add(target1);
				allTargets.add(target2);
				String targetconsensus = nm.getConsensusName(allTargets);
				
				// non-void overlapping edge
				if (overlap_edge_def.getType() != EdgeOntology.VOID_TYPE && sourceconsensus != null && targetconsensus != null)
				{
					if (!allDiffNodes.containsKey(sourceconsensus))
					{
						allDiffNodes.put(sourceconsensus, new Node(sourceconsensus));
					}
					Node sourceresult = allDiffNodes.get(sourceconsensus);
					
					if (!allDiffNodes.containsKey(targetconsensus))
					{
						allDiffNodes.put(targetconsensus, new Node(targetconsensus));
					}
					Node targetresult = allDiffNodes.get(targetconsensus);
					
					Edge overlapdiff = new Edge(sourceresult, targetresult, overlap_edge_def);
					overlap.addEdge(overlapdiff);
				}
			}
		}
		cleaning.removeRedundantEdges(overlap);
		return overlap;
	}
	

	/**
	 * Return one single node from a collection, assuming that there will only be 1
	 * @param nodes the set of nodes
	 * @return the one node in the set, or an UnsupportedOperationException if there are more than 1
	 */
	private Node getSingleNode(Set<Node> nodes)
	{
		if (nodes == null || nodes.isEmpty())
		{
			return null;
		}
		if (nodes.size() > 1)
		{
			throw new UnsupportedOperationException("This algorithm currently only supports 1-1 node mappings");
		}
		return nodes.iterator().next();
	}

	/**
	 * Return one single edge from a collection, assuming that there will only be 1
	 * @param edges the set of edges
	 * @return the one edge in the set, a void one if the set is empty, or throw an UnsupportedOperationException if there are more than 1
	 */
	private EdgeDefinition getSingleEdge(Set<Edge> edges)
	{
		if (edges.size() > 1)
		{
			throw new UnsupportedOperationException("This algorithm currently only supports 1 edge between two nodes in the original networks");
		}
		if (edges.isEmpty())
		{
			return EdgeDefinition.getVoidEdge();
		}
		return edges.iterator().next();
	}

}
