package be.svlandeg.diffany.core.algorithms;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.networks.Condition;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.EdgeDefinition;
import be.svlandeg.diffany.core.networks.EdgeGenerator;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.networks.OverlappingNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.networks.meta.MetaEdge;
import be.svlandeg.diffany.core.networks.meta.MetaEdgeDefinition;
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
	 * An important parameter is overlapNo_cutoff, which determines the amount of support needed for an edge to be included in the 'consensus' condition network. If it equals the number of condition networks, all networks need to agree on an edge.
	 * However, if it is smaller, e.g. 3 out of 4, there can be one 'outlier' network (potentially a different one for each calculated edge), allowing some noise in the input and creating more robust differential networks.
	 * The overlapNo_cutoff should ideally be somewhere between 50% and 100%, but this choice is determined by the specific use-case / application. Instead of being a percentage, this method requires the support to be expressed
	 * as a minimal number of supporting edges (networks).
	 * 
	 * @param reference the reference network
	 * @param conditionNetworks a set of condition-specific networks (at least 2)
	 * @param eo the tree edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param diffName the name to give to the differential network
	 * @param ID the ID of the resulting network
	 * @param overlapNo_cutoff the minimal number of edges that need to overlap between the condition-specific networks
	 * @param weight_cutoff the minimal weight a differential edge should have to be included
	 * 
	 * @return the differential network between the two
	 */
	protected DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, Set<ConditionNetwork> conditionNetworks, TreeEdgeOntology eo, NodeMapper nm, String diffName, int ID, int overlapNo_cutoff, double weight_cutoff)
	{
		DifferentialNetwork diff = new DifferentialNetwork(diffName, ID, reference, conditionNetworks, nm);
		
		Set<Network> allOriginals = new HashSet<Network>();
		allOriginals.addAll(conditionNetworks);
		allOriginals.add(reference);
		
		List<Set<Node>> allEqualSets = nm.getAllEquals(allOriginals);

		Map<String, Node> allDiffNodes = new HashMap<String, Node>();

		Set<String> roots = eo.retrieveAllSourceRootCats(true);
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
				Map<String, Map<Integer, Set<EdgeDefinition>>> edgesBySemanticRoot = new HashMap<String, Map<Integer, Set<EdgeDefinition>>>();
				for (String root : roots)
				{
					edgesBySemanticRoot.put(root, new HashMap<Integer, Set<EdgeDefinition>>());
				}
				for (Network n : allOriginals)
				{
					Set<Edge> allEdges = n.getAllEdges(example_source, example_target);
					for (String root : roots)
					{
						Set<EdgeDefinition> rootEdges = new HashSet<EdgeDefinition>();
						for (Edge e : allEdges)
						{
							String type = e.getType();
							if (eo.isSourceTypeChildOf(type, root) > -1)
							{
								rootEdges.add(e.getDefinition());
							}
						}
						if (rootEdges.isEmpty())
						{
							rootEdges.add(eg.getVoidEdge(eo.isSymmetricalSourceCat(root)));
						}
						edgesBySemanticRoot.get(root).put(n.getID(), rootEdges);
					}
				}

				for (String root : roots)
				{
					boolean symm = eo.isSymmetricalSourceCat(root);

					List<EdgeDefinition> refEdges = new ArrayList<EdgeDefinition>();
					List<EdgeDefinition> conEdges = new ArrayList<EdgeDefinition>();
					List<Set<Integer>> conSupportingNetworks = new ArrayList<Set<Integer>>();
					
					Map<Integer, Set<EdgeDefinition>> edgesByNetwork = edgesBySemanticRoot.get(root);

					for (int networkID : edgesByNetwork.keySet())
					{
						for (EdgeDefinition e : edgesByNetwork.get(networkID))
						{
							if (networkID != reference.getID() && conEdges.contains(e))
							{
								int position = conEdges.indexOf(e);
								conSupportingNetworks.get(position).add(networkID);
							}
							else if (networkID != reference.getID())
							{
								conEdges.add(e);
								Set<Integer> support = new HashSet<Integer>();
								support.add(networkID);
								conSupportingNetworks.add(support);
							}
							else if (networkID == reference.getID() && ! refEdges.contains(e))
							{
								refEdges.add(e);
							}
						}
					}
					if (conEdges.isEmpty())
					{
						conEdges.add(eg.getVoidEdge(symm));
						conSupportingNetworks.add(new HashSet<Integer>());
					}
					if (refEdges.size() > 1)
					{
						throw new IllegalArgumentException("Found more than 1 reference edge in " + reference.getName() + " for semantic root " + root);
					}
					if (refEdges.isEmpty())
					{
						refEdges.add(eg.getVoidEdge(symm));
					}
					EdgeDefinition refEdge = refEdges.iterator().next();
					
					List<EdgeDefinition> mergedEdges = new ArrayList<EdgeDefinition>();
					mergedEdges.add(refEdge);
					mergedEdges.addAll(conEdges);

					List<EdgeDefinition> cleanedEdges = cleaning.unifyDirection(mergedEdges);
					
					EdgeDefinition cleanedRefEdge = cleanedEdges.get(0);
					
					List<EdgeDefinition> cleanedConEdges = new ArrayList<EdgeDefinition>();
					for (int i = 1; i < cleanedEdges.size(); i++)
					{
						System.out.println(" clean con edge " + cleanedEdges.get(i));
						cleanedConEdges.add(cleanedEdges.get(i));
					}

					// TODO shouldn't we check the consensus != null etc (below) BEFORE calculating all these overlap edges ?!
					System.out.println("calculating " + root + " for " + sourceconsensusName + " - " + targetconsensusName);
					EdgeDefinition diff_edge_def = ec.getDifferentialEdge(cleanedRefEdge, cleanedConEdges, conSupportingNetworks, overlapNo_cutoff, weight_cutoff);


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
	 * The additional parameter refRequired determines whether the reference network should always provide support for the overlap edge, or not.
	 * If not, all networks are treated equal. If true, an overlap edge can never be produced if it does not have support in the reference network.
	 * The reference network is determined by trying to cast the input networks to ReferenceNetwork and selecting the one and only unique result.
	 * 
	 * @param networks a set of networks (at least 2).
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param overlap_name the name to give to the overlapping network
	 * @param ID the ID of the resulting network
	 * @param refRequired whether or not the presence of the edge in the reference network is required for it to be included in the overlap network.
	 * When set to true, this method will raise an IllegalArgumentException when no or more than 1 reference network is found.
	 * @param overlapNo_cutoff the minimal number of edges that need to overlap
	 * @param weight_cutoff the minimal weight that the resulting overlap edges should have to be included
	 * @param minOperator whether or not to take the minimum of the edge weights - if false, the maximum is taken
	 * 
	 * @return the overlapping network between the input networks. 
	 * The network will contain Meta edges that store which input networks provide support, given by the original conditions when the input contains ConditionNetworks, or by artificial Conditions
	 * displaying the name of the Networks.
	 * 
	 * TODO v3.0: expand this algorithm to be able to deal with n-m node mappings
	 */
	protected OverlappingNetwork calculateOverlappingNetwork(Set<Network> networks, TreeEdgeOntology eo, NodeMapper nm, String overlap_name, int ID, int overlapNo_cutoff, boolean refRequired, double weight_cutoff, boolean minOperator)
	{
		OverlappingNetwork overlap = new OverlappingNetwork(overlap_name, ID, networks, nm);

		List<Set<Node>> allEqualSets = nm.getAllEquals(networks);

		Set<Integer> refNetworks = new HashSet<Integer>();
		Map<Integer, Network> allNetworks = new HashMap<Integer, Network>();
		for (Network n : networks)
		{
			allNetworks.put(n.getID(), n);
			if (n instanceof ReferenceNetwork)
			{
				refNetworks.add(n.getID());
			}
		}

		if (refRequired && refNetworks.size() != 1)
		{
			String errormsg = "Please define exactly 1 reference network (" + refNetworks.size() + "found) or change the refRequired parameter to false!";
			throw new IllegalArgumentException(errormsg);
		}

		if (allNetworks.size() == 1)
		{
			Network onlyInput = allNetworks.values().iterator().next();
			overlap.setNodesAndEdges(onlyInput.getNodes(), onlyInput.getEdges());
		}

		Map<String, Node> allDiffNodes = new HashMap<String, Node>();

		Set<String> roots = eo.retrieveAllSourceRootCats(true);
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
				Map<String, Map<Integer, Set<EdgeDefinition>>> edgesBySemanticRoot = new HashMap<String, Map<Integer, Set<EdgeDefinition>>>();
				for (String root : roots)
				{
					edgesBySemanticRoot.put(root, new HashMap<Integer, Set<EdgeDefinition>>());
				}
				for (Network n : networks)
				{
					Set<Edge> allEdges = n.getAllEdges(example_source, example_target);
					for (String root : roots)
					{
						Set<EdgeDefinition> rootEdges = new HashSet<EdgeDefinition>();
						for (Edge e : allEdges)
						{
							String type = e.getType();
							if (eo.isSourceTypeChildOf(type, root) > -1)
							{
								rootEdges.add(e.getDefinition());
							}
						}
						if (!rootEdges.isEmpty())
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

					Map<Integer, Set<EdgeDefinition>> edgesByNetwork = edgesBySemanticRoot.get(root);

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
					Map<EdgeDefinition, Set<Integer>> defs = ec.getOverlapEdge(cleanedEdges, supportingNetworks, overlapNo_cutoff, weight_cutoff, minOperator);

					// non-void overlapping edge
					for (EdgeDefinition def : defs.keySet())
					{
						Set<Integer> supports = defs.get(def);
						if (def.getWeight() > 0 && sourceconsensusID != null && targetconsensusID != null)
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

							Set<Condition> conditions = new HashSet<Condition>();
							boolean inReference = false;
							int support = 0;

							if (supports.isEmpty())
							{
								System.out.println("Found no support for " + sourceresult + " - " + targetresult + " - " + def);
							}

							for (int i : supports)
							{
								support++;
								Network input = allNetworks.get(i);
								if (input instanceof ReferenceNetwork)
								{
									inReference = true;
									conditions.add(new Condition(input.getName()));
								}
								else if (input instanceof ConditionNetwork)
								{
									conditions.addAll(((ConditionNetwork) input).getConditions());
								}
								else
								{
									conditions.add(new Condition(input.getName()));
								}
							}
							// only add this edge when the reference network does not matter, or was actually present
							if (!refRequired || inReference)
							{
								MetaEdgeDefinition mergedDef = new MetaEdgeDefinition(def, conditions, support, inReference);

								MetaEdge overlapdiff = new MetaEdge(sourceresult, targetresult, mergedDef);
								overlap.addEdge(overlapdiff);
							}
						}
					}
				}
			}
		}

		return overlap;
	}

}
