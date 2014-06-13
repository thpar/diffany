package be.svlandeg.diffany.core.algorithms;

import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.EdgeDefinition;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.semantics.EdgeOntology;
import be.svlandeg.diffany.core.semantics.NodeMapper;

/**
 * This class provides generic methods useful for network cleaning before or
 * after applying differential algorithms.
 * 
 * @author Sofie Van Landeghem
 */
public class NetworkCleaning
{

	private Logger logger;

	/**
	 * Create a new cleaning object, which can log important messages.
	 * 
	 * @param logger the logger object
	 */
	public NetworkCleaning(Logger logger)
	{
		if (logger == null)
		{
			String errormsg = "The logger should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.logger = logger;
	}

	/**
	 * Clean an output network: Remove redundant/duplicate edges in the network.
	 * 
	 * @param net the network that needs cleaning
	 * @param nm the node mapper for the network
	 */
	public void fullOutputCleaning(Network net)
	{
		removeRedundantSymmetricalEdges(net);
	}

	/**
	 * Remove edges in the network that are symmetrical and are represented
	 * twice (source-target and target-source). One of the two is removed only
	 * then when the type, weight and negation are all equal.
	 * 
	 * @param net the network that needs cleaning
	 */
	protected void removeRedundantSymmetricalEdges(Network net)
	{
		logger.log(" Removing redundant symmetrical edges from network " + net.getName());

		Set<Edge> removed_edges = new HashSet<Edge>();

		Set<Edge> oldEdges = new HashSet<Edge>(net.getEdges());

		// remove duplicate symmetrical edges between source-target and target-source
		for (Edge et : oldEdges)
		{
			Node n1 = et.getSource();
			Node n2 = et.getTarget();

			Set<Edge> all_edges = net.getAllEdges(n1, n2);
			for (Edge eb : all_edges)
			{
				// comparing two edges 'et' and 'eb', not removed in a previous iteration
				if (!et.equals(eb) && !removed_edges.contains(et) && !removed_edges.contains(eb))
				{
					// both are symmetrical and have the same type
					if ((et.isSymmetrical() && eb.isSymmetrical()) && (et.getType().equals(eb.getType())))
					{
						// both have the same weight and negation status
						if ((et.getWeight() == eb.getWeight()) && (et.isNegated() == eb.isNegated()))
						{
							// remove one of the two
							net.removeEdge(eb);
							removed_edges.add(eb);
						}
					}
				}
			}
		}
	}

	/**
	 * Create a new Set of edges with all interactions directed.
	 * 
	 * @param oldSet the old edgeset which might have a mixture of directed and symmetrical edges
	 * 
	 * @return the new edge set with only directed, or only symmetrical edges
	 */
	public Set<EdgeDefinition> makeAllDirected(Set<EdgeDefinition> oldSet)
	{
		Set<EdgeDefinition> newSet = new HashSet<EdgeDefinition>();

		for (EdgeDefinition referenceEdge : oldSet)
		{
			EdgeDefinition newEdge = new EdgeDefinition(referenceEdge);
			newEdge.makeSymmetrical(false);
			newSet.add(newEdge);
		}
		return newSet;
	}

	/**
	 * Clean an input condition-specific network:
	 * Per pair of nodes, group all input edges into subclasses per root category of the EdgeOntology, unify the directionality
	 * (either all symmetric or all directed, as dicated by the edge ontology), and resolve conflicts within a root category.
	 * 
	 * @param net the network that needs cleaning
	 * @param nm the node mapper
	 * @param eo the edge ontology
	 * 
	 * @return a cleaned condition-specific network representing the same semantic information
	 */
	public ConditionNetwork fullInputConditionCleaning(ConditionNetwork net, NodeMapper nm, EdgeOntology eo)
	{
		ConditionNetwork resultNet = new ConditionNetwork(net.getName(), net.getConditions(), nm);
		resultNet.setNodesAndEdges(net.getNodes(), net.getEdges());
		fullCleaning(resultNet, nm, eo);

		return resultNet;
	}

	/**
	 * Clean an input reference network:
	 * Per pair of nodes, group all input edges into subclasses per root category of the EdgeOntology, unify the directionality
	 * (either all symmetric or all directed, as dicated by the edge ontology), and resolve conflicts within a root category.
	 * 
	 * @param net the network that needs cleaning
	 * @param nm the node mapper
	 * @param eo the edge ontology
	 * 
	 * @return a cleaned reference network representing the same semantic information
	 */
	public ReferenceNetwork fullInputRefCleaning(ReferenceNetwork net, NodeMapper nm, EdgeOntology eo)
	{
		ReferenceNetwork resultNet = new ReferenceNetwork(net.getName(), nm);
		resultNet.setNodesAndEdges(net.getNodes(), net.getEdges());
		fullCleaning(resultNet, nm, eo);

		return resultNet;
	}

	/**
	 * Clean a generic input network:
	 * Per pair of nodes, group all input edges into subclasses per root category of the EdgeOntology, unify the directionality
	 * (either all symmetric or all directed, as dictated by the edge ontology), and resolve conflicts within a root category.
	 * 
	 * @param net the network that needs cleaning
	 * @param nm the node mapper
	 * @param eo the edge ontology
	 * 
	 * @return a cleaned reference network representing the same semantic information
	 */
	public InputNetwork fullInputCleaning(InputNetwork net, NodeMapper nm, EdgeOntology eo)
	{
		InputNetwork resultNet = new InputNetwork(net.getName(), nm);
		resultNet.setNodesAndEdges(net.getNodes(), net.getEdges());
		fullCleaning(resultNet, nm, eo);

		return resultNet;
	}

	/**
	 * Clean an input network:
	 * Per pair of nodes, group all input edges into subclasses per root category of the EdgeOntology, unify the directionality
	 * (either all symmetric or all directed, as dicated by the edge ontology), and resolve conflicts within a root category.
	 * 
	 * Be aware: this function changes the input network object!
	 * 
	 * @param net the network that needs cleaning
	 * @param nm the node mapper
	 * @param eo the edge ontology
	 */
	private void fullCleaning(Network net, NodeMapper nm, EdgeOntology eo)
	{
		// make edges directed when defined as such by the edge ontology
		System.out.println("   Unifying direction " + new Date());
		Set<Node> nodes = net.getNodes();
		Set<Edge> edges = new Unification(logger).unifyEdgeDirection(net.getEdges(), eo);
		net.setNodesAndEdges(nodes, edges);

		// clean edges per semantic category
		System.out.println("   Cleaning edges " + new Date());
		cleanEdges(net, nm, eo);
	}

	/**
	 * For each node pair and each semantic root category, resolve conflicts by calling {@link cleanEdgesBetweenNodes}.
	 * This method also removes redundant symmetrical edges at the end.
	 * 
	 * @param net the network that needs cleaning
	 * @param nm the node mapper
	 * @param eo the edge ontology
	 */
	protected void cleanEdges(Network net, NodeMapper nm, EdgeOntology eo)
	{
		Set<Node> allNodes = net.getNodes();
		Set<Edge> newEdges = new HashSet<Edge>();
		Set<String> roots = eo.retrieveAllSourceRootCats();

		// first, determine all node pairs which are relevant in this network
		Set<Edge> oldEdges = net.getEdges();
		Map<String, Set<String>> pairs = new HashMap<String, Set<String>>();

		Map<String, Node> mappedNodes = new HashMap<String, Node>();

		for (Edge e : oldEdges)
		{
			Node source = e.getSource();
			Node target = e.getTarget();

			String sourceID = source.getID();
			String targetID = target.getID();

			if (!pairs.containsKey(sourceID))
			{
				pairs.put(sourceID, new HashSet<String>());
			}
			pairs.get(sourceID).add(targetID);

			if (!pairs.containsKey(targetID))
			{
				pairs.put(targetID, new HashSet<String>());
			}
			pairs.get(targetID).add(sourceID);

			mappedNodes.put(sourceID, source);
			mappedNodes.put(targetID, target);
		}

		// For each node pair, perform a cleaning step
		for (String ID1 : pairs.keySet())
		{
			Node n1 = mappedNodes.get(ID1);
			Set<String> partners = pairs.get(ID1);

			for (String ID2 : partners)
			{
				Node n2 = mappedNodes.get(ID2);

				Set<EdgeDefinition> edges = net.getAllEdgeDefinitions(n1, n2);
				if (!edges.isEmpty())
				{
					Map<String, EdgeDefinition> cleanedEdges = cleanEdgesBetweenNodes(net, eo, roots, edges, n1, n2);
					for (EdgeDefinition def : cleanedEdges.values())
					{
						newEdges.add(new Edge(n1, n2, def));
					}
				}
			}
		}
		net.setNodesAndEdges(allNodes, newEdges);
		// TODO: the above step introduces redundancy which needs to be cleaned again in the next step ... this should be dealt with more properly!
		System.out.println("   Removing redundant symmetry " + new Date());
		removeRedundantSymmetricalEdges(net);

	}

	/**
	 * Clean an input network:
	 * Group all input edges into subclasses per root category of the EdgeOntology, resolve conflicts within a root category.
	 * 
	 * @param net the network that needs cleaning
	 * @param eo the edge ontology
	 * @param roots the root categories of the edge ontology
	 * @param oldEdges all edges between two nodes, including both directed and symmetrical edges
	 * 
	 * @param source the source node (used for logging)
	 * @param target the target node (used for logging)
	 * 
	 * @return all edges grouped by semantic root category, with unified directionality, and only 1 edge per network and root cat.
	 */
	protected Map<String, EdgeDefinition> cleanEdgesBetweenNodes(Network net, EdgeOntology eo, Set<String> roots, Set<EdgeDefinition> oldEdges, Node source, Node target)
	{
		Map<String, Set<EdgeDefinition>> mappedNormalEdges = resolveEdgesPerRoot(eo, roots, oldEdges);
		Map<String, EdgeDefinition> mappedSingleEdges = new HashMap<String, EdgeDefinition>();

		for (String rootCat : mappedNormalEdges.keySet())
		{
			Set<EdgeDefinition> normalEdges = mappedNormalEdges.get(rootCat);
			EdgeDefinition singleEdge = resolveToOne(normalEdges, eo, net.getName(), source, target, rootCat);
			mappedSingleEdges.put(rootCat, singleEdge);
		}
		return mappedSingleEdges;
	}

	/**
	 * Group all input edges into subclasses per root category of the EdgeOntology.
	 * 
	 * @param eo the edge ontology
	 * @param roots the root categories of the edge ontology
	 * @param oldEdges the original set of input edges
	 * @return all input edges grouped by edge root category
	 */
	protected Map<String, Set<EdgeDefinition>> resolveEdgesPerRoot(EdgeOntology eo, Set<String> roots, Set<EdgeDefinition> oldEdges)
	{
		Map<String, Set<EdgeDefinition>> mappedEdges = new HashMap<String, Set<EdgeDefinition>>();

		for (EdgeDefinition refE : oldEdges)
		{
			String edgeType = refE.getType();

			int foundRoot = 0;

			for (String root : roots)
			{
				boolean belongsToRoot = eo.isSourceTypeChildOf(edgeType, root) >= 0;
				if (belongsToRoot)
				{
					foundRoot++;
					if (!mappedEdges.containsKey(root))
					{
						mappedEdges.put(root, new HashSet<EdgeDefinition>());
					}
					mappedEdges.get(root).add(refE);
				}
			}
			if (foundRoot == 0)
			{
				String errorMsg = " Edge source type '" + edgeType + "' could not be linked to a semantic root category in the edge ontology";
				throw new IllegalArgumentException(errorMsg);
			}
			if (foundRoot > 1)
			{
				String errorMsg = " Edge source type '" + edgeType + "' could be linked to more than one semantic root category in the edge ontology";
				throw new IllegalArgumentException(errorMsg);
			}
		}

		return mappedEdges;
	}

	/**
	 * Resolve a set of edges to one. This is currently implemented by taking the edge with the highest weight.
	 * 
	 * It is assumed that resolveEdgesPerRoot was previously used to provide an set of edges which only contains edges for one root category,
	 * and that all edges within this category are either symmetrical, or all directed.
	 * 
	 * @param edges the original set of input edges
	 * @param eo the edge ontology
	 * @param network_name the name/type of input network (used for logging)
	 * @param source the source node (used for logging)
	 * @param target the target node (used for logging)
	 * @param rootCat the root category (used for logging)
	 * 
	 * @return one edge, produced after resolving conflicts, or throws a RunTimeException if no best edge could be found.
	 */
	protected EdgeDefinition resolveToOne(Set<EdgeDefinition> edges, EdgeOntology eo, String network_name, Node source, Node target, String rootCat)
	{
		// TODO v2.0: should we also take into account whether or not one of the edges is more specific / negated?
		double maxWeight = 0.0;
		int numberOriginal = 0;

		for (EdgeDefinition e : edges)
		{
			double weight = e.getWeight();
			if (weight > 0)
			{
				numberOriginal++;
			}
			maxWeight = Math.max(maxWeight, weight);
		}

		for (EdgeDefinition e : edges)
		{
			if (maxWeight == e.getWeight())
			{
				if (numberOriginal > 1)
				{
					logger.log("  Selected only the edge with the highest weight (" + maxWeight + ") between " + source + " and " + target + " for the category " + rootCat + " in " + network_name);
				}
				return e;
			}
		}
		String msg = "Could not resolve the set of edges to one.";
		logger.log("Fatal error: " + msg);
		throw new RuntimeException(msg);
	}

}
