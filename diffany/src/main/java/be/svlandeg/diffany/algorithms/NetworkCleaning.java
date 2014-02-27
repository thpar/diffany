package be.svlandeg.diffany.algorithms;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.EdgeDefinition;
import be.svlandeg.diffany.concepts.Logger;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Node;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
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
		if (log == null)
		{
			String errormsg = "The logger should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.log = log;
	}

	/**
	 * Clean an output network: Remove redundant/duplicate edges in the network.
	 * 
	 * @param net the network that needs cleaning
	 * @param nm the node mapper for the network
	 * @param toLog whether or not to log important messages
	 */
	public void fullOutputCleaning(Network net, boolean toLog)
	{
		removeRedundantSymmetricalEdges(net, toLog);
	}

	/**
	 * Remove edges in the network that are symmetrical and are represented
	 * twice (source-target and target-source). One of the two is removed only
	 * then when the type, weight and negation are all equal.
	 * 
	 * @param net the network that needs cleaning
	 * @param toLog whether or not to log important messages
	 */
	protected void removeRedundantSymmetricalEdges(Network net, boolean toLog)
	{
		if (toLog)
		{
			log.log(" Removing redundant symmetrical edges from the output");
		}

		Set<Edge> removed_edges = new HashSet<Edge>();

		// remove duplicate symmetrical edges between source-target and target-source
		for (Node n1 : net.getNodes())
		{
			for (Node n2 : net.getNodes())
			{
				Set<Edge> all_edges = net.getAllEdges(n1, n2);
				for (Edge et : all_edges)
				{
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
	 * (either all symmetric or all directed), and resolve conflicts within a root category.
	 * 
	 * @param net the network that needs cleaning
	 * @param nm the node mapper
	 * @param eo the edge ontology
	 * @param toLog whether or not to log important messages
	 * 
	 * @return a cleaned condition-specific network representing the same semantic information
	 */
	public ConditionNetwork fullInputConditionCleaning(ConditionNetwork net, NodeMapper nm, EdgeOntology eo, boolean toLog)
	{
		ConditionNetwork resultNet = new ConditionNetwork(net.getName(), net.getConditions(), nm);
		Set<Edge> edges = cleanEdges(net, nm, eo, toLog);
		Set<Node> nodes = net.getNodes();
		resultNet.setNodesAndEdges(nodes, edges);
		return resultNet;
	}
	
	/**
	 * Clean an input reference network:
	 * Per pair of nodes, group all input edges into subclasses per root category of the EdgeOntology, unify the directionality 
	 * (either all symmetric or all directed), and resolve conflicts within a root category.
	 * 
	 * @param net the network that needs cleaning
	 * @param nm the node mapper
	 * @param eo the edge ontology
	 * @param toLog whether or not to log important messages
	 * 
	 * @return a cleaned reference network representing the same semantic information
	 */
	public ReferenceNetwork fullInputRefCleaning(ReferenceNetwork net, NodeMapper nm, EdgeOntology eo, boolean toLog)
	{
		ReferenceNetwork resultNet = new ReferenceNetwork(net.getName(), nm);
		Set<Edge> edges = cleanEdges(net, nm, eo, toLog);
		Set<Node> nodes = net.getNodes();
		resultNet.setNodesAndEdges(nodes, edges);
		return resultNet;
	}
	
	/**
	 * 
	 * @param net
	 * @param nm
	 * @param eo
	 * @param toLog
	 * @return
	 */
	protected Set<Edge> cleanEdges(Network net, NodeMapper nm, EdgeOntology eo, boolean toLog)
	{
		Set<Node> allNodes = net.getNodes();
		Set<Edge> resultEdges = new HashSet<Edge>();
		for (Node source : allNodes)
		{
			for (Node target : allNodes)
			{
				Set<EdgeDefinition> edges = net.getAllEdgeDefinitions(source, target);
				Map<String, EdgeDefinition> cleanedEdges = cleanEdgesBetweenNodes(net, eo, edges, source, target, toLog);
				for (EdgeDefinition def : cleanedEdges.values())
				{
					resultEdges.add(new Edge(source, target, def));
				}
			}
		}
		return resultEdges;
	}
	
	

	/**
	 * Clean an input network:
	 * Group all input edges into subclasses per root category of the EdgeOntology, unify the directionality 
	 * (either all symmetric or all directed), and resolve conflicts within a root category.
	 * 
	 * @param net the network that needs cleaning
	 * @param eo the edge ontology
	 * @param oldEdges all edges between two nodes, including both directed and symmetrical edges
	 * 
	 * @param source the source node (used for logging)
	 * @param target the target node (used for logging)
	 * @param toLog whether or not to log important messages
	 * 
	 * @return all edges grouped by semantic root category, with unified directionality, and only 1 edge per network and root cat.
	 */
	public Map<String, EdgeDefinition> cleanEdgesBetweenNodes(Network net, EdgeOntology eo, Set<EdgeDefinition> oldEdges, Node source, Node target, boolean toLog)
	{
		Map<String, Set<EdgeDefinition>> mappedNormalEdges = resolveEdgesPerRoot(eo, oldEdges);
		Map<String, EdgeDefinition> mappedSingleEdges = new HashMap<String, EdgeDefinition>();

		for (String rootCat : mappedNormalEdges.keySet())
		{
			Set<EdgeDefinition> normalEdges = mappedNormalEdges.get(rootCat);
			EdgeDefinition singleEdge = resolveToOne(normalEdges, eo, net.getName(), source, target, rootCat, toLog);
			mappedSingleEdges.put(rootCat, singleEdge);
		}
		return mappedSingleEdges;
	}

	/**
	 * Group all input edges into subclasses per root category of the EdgeOntology.
	 *
	 * @param eo the edge ontology
	 * @param oldEdges the original set of input edges
	 * @return all input edges grouped by edge root category 
	 */
	protected Map<String, Set<EdgeDefinition>> resolveEdgesPerRoot(EdgeOntology eo, Set<EdgeDefinition> oldEdges)
	{
		Map<String, Set<EdgeDefinition>> mappedEdges = new HashMap<String, Set<EdgeDefinition>>();
		Set<String> roots = eo.retrieveAllSourceRootCats();

		for (EdgeDefinition refE : oldEdges)
		{
			String edgeType = refE.getType();
			String edgeClass = eo.getSourceCategory(edgeType);
			
			int foundRoot = 0;
			
			for (String root : roots)
			{
				boolean belongsToRoot = eo.isSourceChildOf(edgeClass, root) >= 0;
				if (belongsToRoot)
				{
					foundRoot++;
					if (! mappedEdges.containsKey(root))
					{
						mappedEdges.put(root, new HashSet<EdgeDefinition>());
					}
					mappedEdges.get(root).add(refE);
				}
			}
			if (foundRoot == 0)
			{
				String errorMsg = " Edge source type " + edgeType + " could not be linked to a semantic root category in the edge ontology";
				throw new IllegalArgumentException(errorMsg);
			}
			if (foundRoot > 1)
			{
				String errorMsg = " Edge source type " + edgeType + " could be linked to more than one semantic root category in the edge ontology";
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
	 * @param toLog whether or not to log important messages
	 * 
	 * @return one edge, produced after resolving conflicts, or throws a RunTimeException if no best edge could be found. 
	 */
	protected EdgeDefinition resolveToOne(Set<EdgeDefinition> edges, EdgeOntology eo, String network_name, Node source, Node target, String rootCat, boolean toLog)
	{
		// TODO v2.0: should we also take into account whether or not one of the edges is more specific?
		double maxWeight = 0.0;

		for (EdgeDefinition e : edges)
		{
			maxWeight = Math.max(maxWeight, e.getWeight());
		}

		for (EdgeDefinition e : edges)
		{
			if (maxWeight == e.getWeight())
			{
				if (toLog)
				{
					log.log(" Selected only the edge with the highest weight (" + maxWeight + ") between " + source.getName() + " and " + target.getName() + " for the category " + rootCat + " in " + network_name);
				}
				return e;
			}
		}
		String msg = "Could not resolve the set of edges to one.";
		log.log("Fatal error: " + msg);
		throw new RuntimeException(msg);
	}

}
