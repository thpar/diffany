package be.svlandeg.diffany.algorithms;

import java.util.*;

import be.svlandeg.diffany.concepts.*;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * This class can calculate an overlap network between two networks,
 * and is being used by CalculateDiffOfMore for calculating overlap between more than 2 networks
 * 
 * Currently this algorithm assumes a 1-to-1 mapping of nodes between the two networks!
 * 
 * @author Sofie Van Landeghem
 */
public class CalculateDiffOfTwo
{
	
	protected NetworkCleaning cleaning;
	protected Logger log;
	
	/**
	 * Constructor, which initializes the functionality for cleaning a network.
	 * @param log the logger that records logging messages
	 */
	public CalculateDiffOfTwo(Logger log)
	{
		cleaning = new NetworkCleaning(log);
		this.log = log;
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
				Set<Edge> referenceEdges = n1.getAllEdges(source1, target1, nm);
				Set<Edge> back_referenceEdges = n1.getDirectedEdges(target1, source1, nm);

				// get the condition-specific edge
				Set<Edge> conditionEdges = new HashSet<Edge>();
				Set<Edge> back_conditionEdges = new HashSet<Edge>();
				if (source2 != null && target2 != null)
				{
					conditionEdges = n2.getAllEdges(source2, target2, nm);
					back_conditionEdges = n2.getDirectedEdges(target2, source2, nm);
				}
				
				ArrayList<Set<Edge>> condlist = new ArrayList<Set<Edge>>();
				ArrayList<Set<Edge>> back_condlist = new ArrayList<Set<Edge>>();
				condlist.add(conditionEdges);
				back_condlist.add(back_conditionEdges);
				EdgeSet es = new EdgeSet(referenceEdges, condlist);
				EdgeSet back_es = new EdgeSet(back_referenceEdges, back_condlist);
				Map<String, SingleEdgeSet> edgeSets = cleaning.fullInputCleaning(eo, es, back_es);
				
				for (String root : edgeSets.keySet())
				{
					SingleEdgeSet ses = edgeSets.get(root);
					
					Set<EdgeDefinition> allEdges = ses.getAllEdges();
					EdgeDefinition overlap_edge_def = eo.getOverlapEdge(allEdges, cutoff, minOperator);
	
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
		}
		cleaning.fullOutputCleaning(overlap, nm);
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
			String msg = "This algorithm currently only supports 1-1 node mappings.";
			log.log("Fatal error: " + msg);
			throw new UnsupportedOperationException(msg);
		}
		return nodes.iterator().next();
	}

}
