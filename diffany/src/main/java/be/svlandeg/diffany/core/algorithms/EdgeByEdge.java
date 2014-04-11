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
import be.svlandeg.diffany.core.networks.InputNetwork;
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
	
	protected static String EMPTY_NAME = "*empty*";
	
	protected Logger log;
	
	/**
	 * The constructor initializes the algorithm.
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
	 * @param reference the reference network
	 * @param conditions a set of condition-specific networks (at least 2)
	 * @param eo the tree edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param diff_name the name to give to the differential network.
	 * 
	 * @return the differential network between the two
	 */
	protected DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, Set<ConditionNetwork> conditionNetworks, TreeEdgeOntology eo, NodeMapper nm, String diff_name, double cutoff)
	{
		ArrayList<ConditionNetwork> listedConditions = new ArrayList<ConditionNetwork>(conditionNetworks);

		Set<Network> allOriginals = new HashSet<Network>();
		allOriginals.addAll(conditionNetworks);
		allOriginals.add(reference);

		DifferentialNetwork diff = new DifferentialNetwork(diff_name, reference, conditionNetworks, nm);

		Set<Node> allNodes = nm.getAllNodes(allOriginals);

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
					if (!source2.getName().equals(EMPTY_NAME) && !target2.getName().equals(EMPTY_NAME))
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

						String sourceconsensus = nm.getConsensusName(allSources);
						String targetconsensus = nm.getConsensusName(allTargets);

						// non-void differential edge
						if (diff_edge_def.getWeight() > 0 && sourceconsensus != null && targetconsensus != null)
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
			}
		}
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
	 * 
	 * TODO v2.0: calculate overlap directly on set of networks 
	 */
	protected OverlappingNetwork calculateOverlappingNetwork(Set<InputNetwork> networks, TreeEdgeOntology eo, NodeMapper nm, String overlapping_name, double cutoff, boolean minOperator)
	{
		List<InputNetwork> listedNetworks = new ArrayList<InputNetwork>();
		listedNetworks.addAll(networks);

		int numberOfNetworks = listedNetworks.size();
		int first = 0;
		int second = 1;

		Network firstN = listedNetworks.get(first);
		Network secondN = listedNetworks.get(second);
		OverlappingNetwork overlapTmp = calculateOverlappingNetwork(firstN, secondN, eo, nm, overlapping_name, cutoff, minOperator);
		new NetworkCleaning(log).fullOutputCleaning(overlapTmp);
		second++;

		while (second < numberOfNetworks)
		{
			secondN = listedNetworks.get(second);
			overlapTmp = calculateOverlappingNetwork(overlapTmp, secondN, eo, nm, overlapping_name, cutoff, minOperator);
			new NetworkCleaning(log).fullOutputCleaning(overlapTmp);
			second++;
		}
		return overlapTmp;
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
	 * TODO v2.0: calculate overlap directly on set of networks      
	 * TODO v3.0: expand this algorithm to be able to deal with n-m node mappings
	 */
	private OverlappingNetwork calculateOverlappingNetwork(Network n1, Network n2, TreeEdgeOntology eo, NodeMapper nm, String overlap_name, double cutoff, boolean minOperator)
	{
		Set<Network> allOriginals = new HashSet<Network>();
		allOriginals.add(n1);
		allOriginals.add(n2);

		OverlappingNetwork overlap = new OverlappingNetwork(overlap_name, allOriginals, nm);

		Map<Node, Set<Node>> nodeMapping = nm.getAllEquals(n1, n2);
		Set<Node> allNodes = nm.getAllNodes(allOriginals);

		Map<String, Node> allDiffNodes = new HashMap<String, Node>();
		
		Set<String> roots = eo.retrieveAllSourceRootCats();
		EdgeComparison ec = new EdgeComparison(eo);
		EdgeGenerator eg = new EdgeGenerator();

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
				Set<Edge> allRefs = n1.getAllEdges(source1, target1);

				// get the condition-specific edge
				Set<Edge> allConds = new HashSet<Edge>();
				if (source2 != null && target2 != null)
				{
					allConds = n2.getAllEdges(source2, target2);
				}

				for (String root : roots)
				{
					boolean symm = eo.isSymmetricalSourceCat(root);
					Set<EdgeDefinition> rootN1s = new HashSet<EdgeDefinition>();
					for (Edge e : allRefs)
					{
						String type = e.getType();
						if (eo.isSourceTypeChildOf(type, root) > -1)
						{
							rootN1s.add(e);
						}
					}
					if (rootN1s.isEmpty())
					{
						rootN1s.add(eg.getVoidEdge(symm));
					}
					if (rootN1s.size() > 1)
					{
						throw new IllegalArgumentException("Found more than 1 edge in " + n1.getName() + " for semantic root " + root);
					}
					Set<EdgeDefinition> rootN2s = new HashSet<EdgeDefinition>();

					for (Edge e : allConds)
					{
						String type = e.getType();
						if (eo.isSourceTypeChildOf(type, root) > -1)
						{
							rootN2s.add(e);
						}
					}
					if (rootN2s.isEmpty())
					{
						rootN2s.add(eg.getVoidEdge(symm));
					}
					if (rootN2s.size() > 1)
					{
						throw new IllegalArgumentException("Found more than 1 edge in " + n2.getName() + " for semantic root " + root);
					}

					Set<EdgeDefinition> allEdges = new HashSet<EdgeDefinition>();
					allEdges.add(rootN1s.iterator().next());
					allEdges.add(rootN2s.iterator().next());
							
					EdgeDefinition overlap_edge_def = ec.getOverlapEdge(allEdges, cutoff, minOperator);

					Set<Node> allSources = new HashSet<Node>();
					allSources.add(source1);
					allSources.add(source2);
					String sourceconsensus = nm.getConsensusName(allSources);

					Set<Node> allTargets = new HashSet<Node>();
					allTargets.add(target1);
					allTargets.add(target2);
					String targetconsensus = nm.getConsensusName(allTargets);

					// non-void overlapping edge
					if (overlap_edge_def.getWeight() > 0 && sourceconsensus != null && targetconsensus != null)
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
