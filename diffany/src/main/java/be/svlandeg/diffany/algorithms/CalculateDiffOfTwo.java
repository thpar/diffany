package be.svlandeg.diffany.algorithms;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.Node;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * This class can calculate a differential network between one reference and one
 * condition-specific network. 
 * Currently this algorithm assumes a 1-to-1 mapping of nodes between the two networks!
 * 
 * @author Sofie Van Landeghem
 */
public class CalculateDiffOfTwo
{

	/**
	 * Calculate the differential network between the reference and
	 * condition-specific network.
	 * 
	 * @param reference
	 *            the reference network
	 * @param condition
	 *            a condition-specific network
	 * @param eo
	 *            the edge ontology that provides meaning to the edge types
	 * @param nm
	 *            the node mapper that allows to map nodes from the one network
	 *            to the other
	 * @param diff_name
	 *            the name to give to the differential network
	 * @return the differential network between the two
	 * 
	 * TODO: expand this algorithm to be able to deal with n-m node mappings
	 * TODO: expand this algorithm to be able to deal with more than 1 edge between two nodes in the original networks
	 */
	public DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, ConditionNetwork condition, EdgeOntology eo, NodeMapper nm, String diff_name)
	{
		Set<ConditionNetwork> conditions = new HashSet<ConditionNetwork>();
		DifferentialNetwork diff = new DifferentialNetwork(diff_name, reference, conditions);
		Map<Node, Set<Node>> nodeMapping = nm.getAllEquals(reference, condition);
		
		for (Edge edge1 : reference.getEdges())
		{
			Node source1 = edge1.getSource();
			Node target1 = edge1.getTarget();
			
			Set<Node> sources2 = nodeMapping.get(source1);
			if (sources2.size() > 1)
			{
				throw new UnsupportedOperationException("This algorithm currently only supports 1-1 node mappings");
			}
			Node source2 = sources2.iterator().next();
			
			Set<Node> targets2 = nodeMapping.get(target1);
			if (targets2.size() > 1)
			{
				throw new UnsupportedOperationException("This algorithm currently only supports 1-1 node mappings");
			}
			Node target2 = targets2.iterator().next();
			
			Set<Edge> conditionEdges = condition.getAllEdges(source2, target2);
			if (conditionEdges.size() > 1)
			{
				throw new UnsupportedOperationException("This algorithm currently only supports 1 edge between two nodes in the original networks");
			}
			Edge edge2 = conditionEdges.iterator().next();
			
			String diff_edge_category = eo.getDifferentialCategory(edge1, edge2);
			if (diff_edge_category != null) 	// if it's null, there should not be an edge here
			{
				Node sourcediff = new Node(nm.getConsensusName(source1, source2));
				Node targetdiff = new Node(nm.getConsensusName(target1, target2));
				boolean symmetrical = eo.isSymmetrical(diff_edge_category);
				Edge edgediff = new Edge(diff_edge_category, sourcediff, targetdiff, symmetrical);
				diff.addEdge(edgediff);
			}
			
		}
		
		return diff;
	}
}
