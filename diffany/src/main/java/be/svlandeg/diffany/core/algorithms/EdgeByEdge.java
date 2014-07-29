package be.svlandeg.diffany.core.algorithms;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.EdgeDefinition;
import be.svlandeg.diffany.core.networks.EdgeGenerator;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.networks.OverlappingNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

/**
 * This class calculates overlap/differential networks on an edge-by-edge basis.
 * 
 * @author Sofie Van Landeghem
 */
public class EdgeByEdge
{

	protected static String EMPTY_ID = "*empty*";
	protected static String EMPTY_DISPLAY_NAME = "Empty node";

	protected Logger log;

	/**
	 * The constructor initializes the algorithm.
	 * 
	 * @param log the logger that records logging messages
	 */
	public EdgeByEdge(Logger log)
	{
		this.log = log;
	}

	/**
	 * Calculate the differential network between the reference and condition-specific networks.
	 * The overlapping network should be calculated independently!
	 * This method can only be called from within the package (CalculateDiff) and can thus assume proper input.
	 * 
	 * TODO: currently the differential network is calculated with the philosophy of 100% overlap between the condition-specific networks. Can this definition be made more fuzzy?
	 * 
	 * @param reference the reference network
	 * @param conditions a set of condition-specific networks (at least 2)
	 * @param eo the tree edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param diff_name the name to give to the differential network
	 * @param ID the ID of the resulting network
	 * 
	 * @return the differential network between the two
	 */
	protected DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, Set<ConditionNetwork> conditionNetworks, TreeEdgeOntology eo, NodeMapper nm, String diff_name, int ID, double cutoff)
	{
		ArrayList<ConditionNetwork> listedConditions = new ArrayList<ConditionNetwork>(conditionNetworks);

		Set<Network> allOriginals = new HashSet<Network>();
		allOriginals.addAll(conditionNetworks);
		allOriginals.add(reference);

		DifferentialNetwork diff = new DifferentialNetwork(diff_name, ID, reference, conditionNetworks, nm);

		Set<Node> allNodes = nm.getAllNodes(allOriginals);

		// nodes in the differential network, by ID
		Map<String, Node> allDiffNodes = new HashMap<String, Node>();

		Set<String> roots = eo.retrieveAllSourceRootCats();
		EdgeComparison ec = new EdgeComparison(eo);
		EdgeGenerator eg = new EdgeGenerator();

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
					allsources2.add(new Node(EMPTY_ID, EMPTY_DISPLAY_NAME));
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
						alltargets2.add(new Node(EMPTY_ID, EMPTY_DISPLAY_NAME));
					}
				}

				// get the reference edges
				Set<Edge> allRefs = reference.getAllEdges(source1, target1);

				// get all condition-specific edges
				ArrayList<Set<Edge>> condlist = new ArrayList<Set<Edge>>();

				for (int i = 0; i < listedConditions.size(); i++)
				{
					Set<Edge> conditionEdges = new HashSet<Edge>();

					ConditionNetwork condition = listedConditions.get(i);
					Node source2 = allsources2.get(i);
					Node target2 = alltargets2.get(i);
					if (!source2.getID().equals(EMPTY_ID) && !target2.getID().equals(EMPTY_ID))
					{
						conditionEdges = condition.getAllEdges(source2, target2);
					}
					condlist.add(i, conditionEdges);
				}

				for (String root : roots)
				{
					boolean aRef = true;
					boolean symm = eo.isSymmetricalSourceCat(root);
					Set<EdgeDefinition> rootRefs = new HashSet<EdgeDefinition>();
					for (Edge e : allRefs)
					{
						String type = e.getType();
						if (eo.isSourceTypeChildOf(type, root) > -1)
						{
							rootRefs.add(e);
						}
					}
					if (rootRefs.isEmpty())
					{
						aRef = false;
						rootRefs.add(eg.getVoidEdge(symm));
					}
					if (rootRefs.size() > 1)
					{
						throw new IllegalArgumentException("Found more than 1 reference edge in " + reference.getName() + " for semantic root " + root);
					}
					Set<EdgeDefinition> rootCons = new HashSet<EdgeDefinition>();
					boolean atLeastOneCon = false;
					for (int i = 0; i < listedConditions.size(); i++)
					{
						Set<EdgeDefinition> thisRootCons = new HashSet<EdgeDefinition>();
						for (Edge e : condlist.get(i))
						{
							String type = e.getType();
							if (eo.isSourceTypeChildOf(type, root) > -1)
							{
								thisRootCons.add(e);
							}
						}
						if (thisRootCons.isEmpty())
						{
							thisRootCons.add(eg.getVoidEdge(symm));
						}
						else if (thisRootCons.size() > 1)
						{
							throw new IllegalArgumentException("Found more than 1 condition edge in " + listedConditions.get(i).getName() + " for semantic root " + root);

						}
						else
						{
							atLeastOneCon = true;
						}
						rootCons.add(thisRootCons.iterator().next());
					}
					if (atLeastOneCon || aRef)
					{
						EdgeDefinition diff_edge_def = ec.getDifferentialEdge(rootRefs.iterator().next(), rootCons, cutoff);

						String sourceconsensusName = nm.getConsensusName(allSources);
						String sourceconsensusID = nm.getConsensusID(allSources);
						String targetconsensusName = nm.getConsensusName(allTargets);
						String targetconsensusID = nm.getConsensusID(allTargets);

						// non-void differential edge
						if (diff_edge_def.getWeight() > 0 && sourceconsensusID != null && targetconsensusID != null)
						{
							if (!allDiffNodes.containsKey(sourceconsensusID))
							{
								allDiffNodes.put(sourceconsensusID, new Node(sourceconsensusID, sourceconsensusName));
							}
							Node sourceresult = allDiffNodes.get(sourceconsensusID);

							if (!allDiffNodes.containsKey(targetconsensusID))
							{
								allDiffNodes.put(targetconsensusID, new Node(targetconsensusID, targetconsensusName));
							}
							Node targetresult = allDiffNodes.get(targetconsensusID);

							Edge edgediff = new Edge(sourceresult, targetresult, diff_edge_def);
							diff.addEdge(edgediff);
						}
					}
				}
			}
		}
		return diff;
	}

	/**
	 * Calculate the overlapping network between a set of networks.
	 * This method can only be called from within the package (CalculateDiff) and can thus assume proper input.
	 * 
	 * An important parameter is overlapNo_cutoff, which determines the amount of support needed for an edge to be included in the overlap network. If it equals the number of input networks, all networks need to agree on an edge.
	 * However, if it is smaller, e.g. 3 out of 4, there can be one 'outlier' network (potentially a different one for each calculated edge), allowing some noise in the input and creating more robust overlap networks.
	 * The overlapNo_cutoff should ideally be somewhere between 50% and 100%, but this choice is determined by the specific use-case / application. Instead of being a percentage, this method requires the support to be expressed
	 * as a minimal number of supporting edges (networks).
	 * 
	 * @param networks a set of networks (at least 2)
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param overlap_name the name to give to the overlapping network
	 * @param ID the ID of the resulting network
	 * @param overlapNo_cutoff the minimal number of edges that need to overlap
	 * @param weight_cutoff the minimal weight that the resulting overlap edges should have to be included
	 * @param minOperator whether or not to take the minimum of the edge weights - if false, the maximum is taken
	 * 
	 * @return the overlapping network between the two
	 * 
	 * TODO v3.0: expand this algorithm to be able to deal with n-m node mappings
	 */
	protected OverlappingNetwork calculateOverlappingNetwork(Set<Network> networks, TreeEdgeOntology eo, NodeMapper nm, String overlap_name, int ID, int overlapNo_cutoff, double weight_cutoff, boolean minOperator)
	{
		List<Network> listedNetworks = new ArrayList<Network>();
		listedNetworks.addAll(networks);

		OverlappingNetwork overlap = new OverlappingNetwork(overlap_name, ID, networks, nm);

		List<Set<Node>> allEqualSets = nm.getAllEquals(networks);

		Map<String, Node> allDiffNodes = new HashMap<String, Node>();

		Set<String> roots = eo.retrieveAllSourceRootCats();
		EdgeComparison ec = new EdgeComparison(eo);
		EdgeGenerator eg = new EdgeGenerator();
		NetworkCleaning cleaning = new NetworkCleaning(log);

		for (Set<Node> sources : allEqualSets) // source nodes (equals across networks) 
		{
			Node example_source = sources.iterator().next();
			String sourceconsensusID = nm.getConsensusID(sources);
			String sourceconsensusName = nm.getConsensusName(sources);

			for (Set<Node> targets : allEqualSets) // target nodes (equals across networks)
			{
				Node example_target = targets.iterator().next();
				String targetconsensusID = nm.getConsensusID(targets);
				String targetconsensusName = nm.getConsensusName(targets);

				// get all edges from the input network
				Map<String, Map<Integer, Set<Edge>>> edgesBySemanticRoot = new HashMap<String, Map<Integer, Set<Edge>>>();
				for (String root : roots)
				{
					edgesBySemanticRoot.put(root, new HashMap<Integer, Set<Edge>>());
				}
				for (Network n : networks)
				{
					Set<Edge> allEdges = n.getAllEdges(example_source, example_target);
					for (String root : roots)
					{
						Set<Edge> rootEdges = new HashSet<Edge>();
						for (Edge e : allEdges)
						{
							String type = e.getType();
							if (eo.isSourceTypeChildOf(type, root) > -1)
							{
								rootEdges.add(e);
							}
						}
						if (! rootEdges.isEmpty())
						{
							edgesBySemanticRoot.get(root).put(n.getID(), rootEdges);
						}
					}
					
				}

				for (String root : roots)
				{
					boolean symm = eo.isSymmetricalSourceCat(root);
					
					List<EdgeDefinition> edges = new ArrayList<EdgeDefinition>();
					List<Set<Integer>> supportingNetworks = new ArrayList<Set<Integer>>();
					
					Map<Integer, Set<Edge>> edgesByNetwork = edgesBySemanticRoot.get(root);
					
					for (int networkID : edgesByNetwork.keySet())
					{
						for (EdgeDefinition e : edgesByNetwork.get(networkID))
						{
							if (edges.contains(e))
							{
								int position = edges.indexOf(e);
								supportingNetworks.get(position).add(networkID);
							}
							else
							{
								edges.add(e);
								Set<Integer> support = new HashSet<Integer>();
								support.add(networkID);
								supportingNetworks.add(support);
							}
						}
					}

					if (edges.isEmpty())
					{
						edges.add(eg.getVoidEdge(symm));
						supportingNetworks.add(new HashSet<Integer>());
					}

					List<EdgeDefinition> cleanedEdges = cleaning.unifyDirection(edges);
					
					// TODO shouldn't we check the consensus != null etc (below) BEFORE calculating all these overlap edges ?!
					Map<EdgeDefinition, Set<Integer>> overlap_edge_defs = ec.getOverlapEdge(cleanedEdges, supportingNetworks, overlapNo_cutoff, weight_cutoff, minOperator);

					// non-void overlapping edge
					for (EdgeDefinition overlap_edge_def : overlap_edge_defs.keySet())
					{
						if (overlap_edge_def.getWeight() > 0 && sourceconsensusID != null && targetconsensusID != null)
						{
							if (!allDiffNodes.containsKey(sourceconsensusID))
							{
								allDiffNodes.put(sourceconsensusID, new Node(sourceconsensusID, sourceconsensusName));
							}
							Node sourceresult = allDiffNodes.get(sourceconsensusID);

							if (!allDiffNodes.containsKey(targetconsensusID))
							{
								allDiffNodes.put(targetconsensusID, new Node(targetconsensusID, targetconsensusName));
							}
							Node targetresult = allDiffNodes.get(targetconsensusID);

							Edge overlapdiff = new Edge(sourceresult, targetresult, overlap_edge_def);
							overlap.addEdge(overlapdiff);
						}
					}
				}
			}
		}

		return overlap;
	}

	/**
	 * Return one single node from a collection, assuming that there will only be 1
	 * 
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
