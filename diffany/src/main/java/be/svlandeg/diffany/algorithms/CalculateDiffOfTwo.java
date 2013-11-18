package be.svlandeg.diffany.algorithms;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Node;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.concepts.SharedNetwork;
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

	/**
	 * Calculate the differential network between the reference and condition-specific network. 
	 * The differential network will also store theshared/overlap/'house-keeping' interactions.
	 * 
	 * @param reference the reference network
	 * @param condition a condition-specific network
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param diff_name the name to give to the differential network. 
	 * 	The associated shared network will be named similarly with appendix '_overlap'.
	 * @return the differential network between the two
	 *         
	 * TODO: expand this algorithm to be able to deal with n-m node mappings (v.2.0) 
	 * TODO: expand this algorithm to be able to deal with more than 1 edge between two nodes 
	 * in the original networks (v.1.0)
	 * TODO: properly test influence of symmetric edges etc.
	 */
	public DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, ConditionNetwork condition, EdgeOntology eo, NodeMapper nm, String diff_name)
	{
		Set<ConditionNetwork> conditions = new HashSet<ConditionNetwork>();
		conditions.add(condition);

		Set<Network> allOriginals = new HashSet<Network>();
		allOriginals.addAll(conditions);
		allOriginals.add(reference);

		DifferentialNetwork diff = new DifferentialNetwork(diff_name, reference, conditions);
		SharedNetwork shared = new SharedNetwork(diff_name + CalculateDiff.overlapnamesuffix, allOriginals);

		Map<Node, Set<Node>> nodeMapping = nm.getAllEquals(reference, condition);

		Set<Node> allNodes = new HashSet<Node>();
		allNodes.addAll(reference.getNodes());
		allNodes.addAll(condition.getNodes());

		Map<String, Node> allConditionNodes = new HashMap<String, Node>();
		allNodes.addAll(reference.getNodes());
		allNodes.addAll(condition.getNodes());
		
		for (Node source1 : allNodes) // source node in reference network
		{
			// get the equivalent source node in the condition network
			Set<Node> sources2 = nodeMapping.get(source1);
			Node source2 = getSingleNode(sources2);

			for (Node target1 : allNodes) // target node in reference network
			{
				// get the equivalent target node in the condition network
				Set<Node> targets2 = nodeMapping.get(target1);
				Node target2 = getSingleNode(targets2);

				// get the reference edge
				Set<Edge> referenceEdges = reference.getAllEdges(source1, target1, true);
				Edge edge1 = getSingleEdge(referenceEdges);

				// get the condition-specific edge
				Set<Edge> conditionEdges = new HashSet<Edge>();
				if (source2 != null && target2 != null)
				{
					conditionEdges = condition.getAllEdges(source2, target2, true);
				}
				Edge edge2 = getSingleEdge(conditionEdges);

				String diff_edge_category = eo.getDifferentialCategory(getEdgeCat(edge1, eo), getEdgeCat(edge2, eo));

				String sourceconsensus = nm.getConsensusName(source1, source2);
				if (! allConditionNodes.containsKey(sourceconsensus))
				{
					allConditionNodes.put(sourceconsensus, new Node(sourceconsensus));
				}
				Node sourceresult = allConditionNodes.get(sourceconsensus);
				
				String targetconsensus = nm.getConsensusName(target1, target2);
				if (! allConditionNodes.containsKey(targetconsensus))
				{
					allConditionNodes.put(targetconsensus, new Node(targetconsensus));
				}
				Node targetresult = allConditionNodes.get(targetconsensus);
				
				// overlapping edge: either same category, or both void
				if (diff_edge_category == EdgeOntology.VOID_EDGE)
				{
					String edge_cat = getEdgeCat(edge1, eo);

					// only add to the network if there was in fact an edge in the reference network
					if (edge_cat != EdgeOntology.VOID_EDGE)
					{
						boolean symmetrical = edge1.isSymmetrical();
						Edge edgeoverlap = new Edge(edge_cat, sourceresult, targetresult, symmetrical);
						shared.addEdge(edgeoverlap);
					}
				}
				else
				{
					boolean symmetrical = eo.isSymmetrical(diff_edge_category);
					Edge edgediff = new Edge(diff_edge_category, sourceresult, targetresult, symmetrical);
					diff.addEdge(edgediff);
				}
			}
		}
		diff.removeRedundantEdges();
		shared.removeRedundantEdges();
		
		diff.setSharedNetwork(shared);
		return diff;
	}

	/**
	 * Return one single node from a collection, assuming that there will only be 1
	 * @param nodes the set of nodes
	 * @return the one node in the set, or an UnsupportedOperationException if there are more than 1
	 */
	private Node getSingleNode(Set<Node> nodes)
	{
		if (nodes.size() > 1)
		{
			throw new UnsupportedOperationException("This algorithm currently only supports 1-1 node mappings");
		}
		if (nodes.isEmpty())
		{
			return null;
		}
		return nodes.iterator().next();
	}

	/**
	 * Return one single edge from a collection, assuming that there will only be 1
	 * @param edges the set of edges
	 * @return the one edge in the set, or throw an UnsupportedOperationException if there are more than 1
	 */
	private Edge getSingleEdge(Set<Edge> edges)
	{
		if (edges.size() > 1)
		{
			throw new UnsupportedOperationException("This algorithm currently only supports 1 edge between two nodes in the original networks");
		}
		if (edges.isEmpty())
		{
			return null;
		}
		return edges.iterator().next();
	}

	/**
	 * Return the edge type of an edge, or EdgeOntology.VOID_EDGE when there is no edge (input parameter is null)
	 * @param edge the edge (may be null)
	 * @return the type of the edge (or EdgeOntology.VOID_EDGE)
	 */
	private String getEdgeCat(Edge edge, EdgeOntology eo)
	{
		if (edge == null)
		{
			return EdgeOntology.VOID_EDGE;
		}
		return eo.getCategory(edge.getType());
	}
}
