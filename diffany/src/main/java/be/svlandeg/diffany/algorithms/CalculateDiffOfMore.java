package be.svlandeg.diffany.algorithms;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.EdgeDefinition;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Node;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.concepts.OverlappingNetwork;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * This class can calculate a differential network between one reference and a set of condition-specific networks. 
 * Currently this algorithm assumes a 1-to-1 mapping of nodes between the two networks!
 * 
 * @author Sofie Van Landeghem
 */
public class CalculateDiffOfMore
{
	
	protected NetworkCleaning cleaning;
	protected CalculateDiffOfTwo twoProcessor;
	protected static String EMPTY_NAME = "*empty*";
	
	/**
	 * Constructor initializes the algorithm suites.
	 */
	public CalculateDiffOfMore()
	{
		twoProcessor = new CalculateDiffOfTwo();
		cleaning = new NetworkCleaning();
	}
	
	/**
	 * Calculate the differential network between the reference and condition-specific networks. 
	 * The overlapping network should be calculated independently!
	 * This method can only be called from within the package (CalculateDiff) and can thus assume proper input.
	 * 
	 * @param reference the reference network
	 * @param conditions a set of condition-specific networks (at least 2)
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param diff_name the name to give to the differential network. 
	 * 
	 * @return the differential network between the two    
	 */
	protected DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, Set<ConditionNetwork> conditionNetworks, EdgeOntology eo, 
			NodeMapper nm, String diff_name, double cutoff)
	{
		ArrayList<ConditionNetwork> listedConditions = new ArrayList<ConditionNetwork>(conditionNetworks);
		
		Set<Network> allOriginals = new HashSet<Network>();
		allOriginals.addAll(conditionNetworks);
		allOriginals.add(reference);

		DifferentialNetwork diff = new DifferentialNetwork(diff_name, reference, conditionNetworks);

		Set<Node> allNodes = nm.getAllNodes(allOriginals);

		Map<String, Node> allDiffNodes = new HashMap<String, Node>();

		for (Node source1 : allNodes) // source node in reference network
		{
			// get the equivalent source nodes in the condition networks
			List<Node> allsources2 = new ArrayList<Node>();
			Set<Node> allSources = new HashSet<Node>();
			allSources.add(source1);
			for (ConditionNetwork condition : listedConditions)
			{
				Map<Node, Set<Node>> nodeMapping = nm.getAllEquals(reference, condition);
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
				if (source2 != null)
				{
					allsources2.add(source2);
					allSources.add(source2);
				}
				else
				{
					allsources2.add(new Node(EMPTY_NAME));
				}
			}
			

			for (Node target1 : allNodes) // target node in reference network
			{
				// get the equivalent target nodes in the condition networks
				List<Node> alltargets2 = new ArrayList<Node>();
				Set<Node> allTargets = new HashSet<Node>();
				allTargets.add(target1);
				for (ConditionNetwork condition : listedConditions)
				{
					Map<Node, Set<Node>> nodeMapping = nm.getAllEquals(reference, condition);
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
					if (target2 != null)
					{
						alltargets2.add(target2);
						allTargets.add(target2);
					}
					else
					{
						alltargets2.add(new Node(EMPTY_NAME));
					}
				}

				// get the reference edge
				Set<Edge> referenceEdges = reference.getAllEdges(source1, target1);
				EdgeDefinition edgedef1 = getSingleEdge(referenceEdges);

				// get all condition-specific edges (one for each condition network)
				Set<EdgeDefinition> edgesdefs2 = new HashSet<EdgeDefinition>();
				for (int i = 0; i < listedConditions.size(); i++)
				{
					ConditionNetwork condition = listedConditions.get(i);
					Node source2 = allsources2.get(i);
					Node target2 = alltargets2.get(i);
					if (! source2.getName().equals(EMPTY_NAME) && ! target2.getName().equals(EMPTY_NAME))
					{
						Set<Edge> edges2 = condition.getAllEdges(source2, target2);
						edgesdefs2.addAll(edges2);
						if (edges2.isEmpty())
						{
							edgesdefs2.add(EdgeDefinition.getVoidEdge());
						}
					}
				}	
				EdgeDefinition diff_edge_def = eo.getDifferentialEdge(edgedef1, edgesdefs2, cutoff);

				String sourceconsensus = nm.getConsensusName(allSources);
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
		cleaning.directSymmetricalWhenOverlapping(diff, eo);
		
		return diff;
	}
	
	
	/**
	 * Calculate the overlapping network between a set of networks. 
	 * This method can only be called from within the package (CalculateDiff) and can thus assume proper input.
	 * 
	 * @param networks a set of networks (at least 2)
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param overlapping_name the name to give to the overlapping network. 
	 * @param minOperator whether or not to take the minimum of the edge weights for the overlapping edges - if false, the maximum is taken
	 * 
	 * @return the differential network between the two
	 */
	protected OverlappingNetwork calculateOverlappingNetwork(Set<Network> networks, EdgeOntology eo, 
			NodeMapper nm, String overlapping_name, double cutoff, boolean minOperator)
	{
		List<Network> listedNetworks = new ArrayList<Network>();
		listedNetworks.addAll(networks);
		
		int numberOfNetworks = listedNetworks.size();
		int first = 0;
		int second = 1;
		
		Network firstN = listedNetworks.get(first);
		Network secondN = listedNetworks.get(second);
		OverlappingNetwork overlapTmp = twoProcessor.calculateOverlappingNetwork(firstN, secondN, eo, nm, overlapping_name, cutoff, minOperator);
		second++;
		
		while (second < numberOfNetworks)
		{	
			secondN = listedNetworks.get(second);
			overlapTmp = twoProcessor.calculateOverlappingNetwork(overlapTmp, secondN, eo, nm, overlapping_name, cutoff, minOperator);
			second++;
		}
		
		return overlapTmp;
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
