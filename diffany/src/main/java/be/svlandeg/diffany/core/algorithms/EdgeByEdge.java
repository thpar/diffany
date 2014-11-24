package be.svlandeg.diffany.core.algorithms;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.networks.Condition;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.EdgeDefinition;
import be.svlandeg.diffany.core.networks.EdgeGenerator;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.progress.ScheduledTask;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

/**
 * This class calculates consensus/differential networks on an edge-by-edge basis.
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
	 * The consensus network should be calculated independently!
	 * This method can only be called from within the package (CalculateDiff) and can thus assume proper input.
	 * 
	 * An important parameter is the supportingCutoff, which determines the amount of support needed for an edge to be included in the 'consensus' condition network. If it equals the number of condition networks, all networks need to agree on an edge.
	 * However, if it is smaller, e.g. 3 out of 4, there can be one 'outlier' network (potentially a different one for each calculated edge), allowing some noise in the input and creating more robust differential networks.
	 * The supportingCutoff should ideally be somewhere between 50% and 100%, but this choice is determined by the specific use-case / application. Instead of being a percentage, this method requires the support to be expressed
	 * as a minimal number of supporting edges (networks).
	 * 
	 * @param reference the reference network
	 * @param conditionNetworks a set of condition-specific networks (at least 2)
	 * @param eo the tree edge ontology that provides meaning to the edge types
	 * @param diffName the name to give to the differential network
	 * @param ID the ID of the resulting network
	 * @param supportingCutoff the minimal number of networks that need to support an edge in the consensus network
	 * @param weightCutoff the minimal weight a differential edge should have to be included
	 * @param task the task object that keeps track of the progress of this calculation (can be null)
	 * 
	 * @return the differential network between the two
	 */
	protected DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, Set<ConditionNetwork> conditionNetworks, TreeEdgeOntology eo, String diffName, int ID, int supportingCutoff, double weightCutoff, ScheduledTask task)
	{
		DifferentialNetwork diff = new DifferentialNetwork(diffName, ID, reference, conditionNetworks);

		Set<Network> allOriginals = new HashSet<Network>();
		allOriginals.add(reference);
		allOriginals.addAll(conditionNetworks);

		Map<String, Set<String>> source2targetIDs = new HashMap<String, Set<String>>();
		Map<String, Set<Node>> condNodesByID = new HashMap<String, Set<Node>>();
		Map<String, Node> refNodesByID = new HashMap<String, Node>();

		for (Network n : conditionNetworks)
		{
			for (Edge e : n.getEdges())
			{
				Node source = e.getSource();
				String sourceID = source.getID();
				Node target = e.getTarget();
				String targetID = target.getID();

				// add the source-target mapping
				if (!source2targetIDs.containsKey(sourceID))
				{
					source2targetIDs.put(sourceID, new HashSet<String>());
				}
				source2targetIDs.get(sourceID).add(targetID);

				// add the source node by ID
				if (!condNodesByID.containsKey(sourceID))
				{
					condNodesByID.put(sourceID, new HashSet<Node>());
				}
				condNodesByID.get(sourceID).add(source);

				// add the target node by ID
				if (!condNodesByID.containsKey(targetID))
				{
					condNodesByID.put(targetID, new HashSet<Node>());
				}
				condNodesByID.get(targetID).add(target);
			}
		}
		
		for (Edge e : reference.getEdges())
		{
			Node source = e.getSource();
			String sourceID = source.getID();
			Node target = e.getTarget();
			String targetID = target.getID();

			// add the source-target mapping
			if (!source2targetIDs.containsKey(sourceID))
			{
				source2targetIDs.put(sourceID, new HashSet<String>());
			}
			source2targetIDs.get(sourceID).add(targetID);

			// add the source node by ID
			refNodesByID.put(sourceID, source);

			// add the target node by ID
			refNodesByID.put(targetID, target);
		}

		Map<String, Node> allDiffNodes = new HashMap<String, Node>();

		Set<String> roots = eo.retrieveAllSourceRootCats(true);

		NodeComparison nc = new NodeComparison();
		EdgeComparison ec = new EdgeComparison(eo);
		EdgeGenerator eg = new EdgeGenerator();
		NetworkCleaning cleaning = new NetworkCleaning(log);

		String conditionNetworkNames = "";
		for (ConditionNetwork c : conditionNetworks)
		{
			conditionNetworkNames += c.getName() + "/";
		}
		conditionNetworkNames = conditionNetworkNames.substring(0, conditionNetworkNames.length() - 1); // cutoff the last "/"

		int totalPairs = 0;
		for (String sourceID : source2targetIDs.keySet())
		{
			totalPairs += source2targetIDs.get(sourceID).size();
		}
		
		int ticksPerReport = 1;
		int pairsPerReport = 0;
		
		if (task != null)
		{
			String progressMessage = "Calculating differential network between " + reference.getName() + " and " + conditionNetworkNames;
			task.setMessage(progressMessage);
			task.ticksDone(ticksPerReport);
			pairsPerReport = ticksPerReport * totalPairs / task.ticksToGo();
		}

		int pairsProgressed = 0;
		for (String sourceID : source2targetIDs.keySet())
		{
			// Determine the node with this source ID in the differential network
			Node diff_source = allDiffNodes.get(sourceID);
			if (diff_source == null)
			{
				diff_source = nc.getConsensusNode(refNodesByID.get(sourceID), condNodesByID.get(sourceID), supportingCutoff);	
				allDiffNodes.put(sourceID, diff_source);
			}

			for (String targetID : source2targetIDs.get(sourceID)) // all possible targets for this source
			{
				pairsProgressed++;
				// Determine the node with this target ID in the differential network
				Node diff_target = allDiffNodes.get(targetID);
				if (diff_target == null)
				{
					diff_target = nc.getConsensusNode(refNodesByID.get(targetID), condNodesByID.get(targetID), supportingCutoff);	
					allDiffNodes.put(targetID, diff_target);
				}

				// notify the progress listener of our progress
				if (task != null && pairsPerReport > 0 && pairsProgressed % pairsPerReport == 0)
				{
					task.ticksDone(ticksPerReport);
				}
				
				// get all edges from the input network
				Map<String, Map<Integer, Set<EdgeDefinition>>> edgesBySemanticRoot = new HashMap<String, Map<Integer, Set<EdgeDefinition>>>();
				for (String root : roots)
				{
					edgesBySemanticRoot.put(root, new HashMap<Integer, Set<EdgeDefinition>>());
				}
				for (Network n : allOriginals)
				{
					Set<Edge> allEdges = n.getAllEdges(sourceID, targetID);
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
							else if (networkID == reference.getID() && !refEdges.contains(e))
							{
								if (! (e.isNegated() || e.getWeight() == 0))
								{
									refEdges.add(e);
								}
							}
						}
					}
					
					/* It only makes sense to try and calculate something if we have at least 1 non-void edge */
					if (!conEdges.isEmpty() || !refEdges.isEmpty())
					{
						if (conEdges.isEmpty())
						{
							conEdges.add(eg.getVoidEdge(symm));
							conSupportingNetworks.add(new HashSet<Integer>());
						}
						if (refEdges.isEmpty())
						{
							refEdges.add(eg.getVoidEdge(symm));
						}
						if (refEdges.size() > 1)
						{
							throw new IllegalArgumentException("Found more than 1 non-zero reference edge in " + reference.getName() + " for semantic root " + root);
						}
						
						List<EdgeDefinition> mergedEdges = new ArrayList<EdgeDefinition>();
						mergedEdges.addAll(refEdges);
						mergedEdges.addAll(conEdges);

						List<EdgeDefinition> cleanedEdges = cleaning.unifyDirection(mergedEdges);

						EdgeDefinition cleanedRefEdge = cleanedEdges.get(0);

						List<EdgeDefinition> cleanedConEdges = new ArrayList<EdgeDefinition>();
						for (int i = 1; i < cleanedEdges.size(); i++)
						{
							cleanedConEdges.add(cleanedEdges.get(i));
						}

						EdgeDefinition diff_edge_def = ec.getDifferentialEdge(cleanedRefEdge, cleanedConEdges, conSupportingNetworks, supportingCutoff, weightCutoff);

						// non-void differential edge
						if (diff_edge_def.getWeight() > 0)
						{
							Edge edgediff = new Edge(diff_source, diff_target, diff_edge_def);
							diff.addEdge(edgediff);
						}
					}
				}
			}
		}
		// notify the progress listener of the fact that we're done (100%)
		if (task != null)
		{
			task.done();
		}

		return diff;
	}
	

	/**
	 * Calculate the consensus network between a set of networks.
	 * This method can only be called from within the package (CalculateDiff) and can thus assume proper input.
	 * 
	 * An important parameter is supportingCutoff, which determines the amount of support needed for an edge to be included in the consensus network. If it equals the number of input networks, all networks need to agree on an edge.
	 * However, if it is smaller, e.g. 3 out of 4, there can be one 'outlier' network (potentially a different one for each calculated edge), allowing some noise in the input and creating more robust consensus networks.
	 * The supportingCutoff should ideally be somewhere between 50% and 100%, but this choice is determined by the specific use-case / application. Instead of being a percentage, this method requires the support to be expressed
	 * as a minimal number of supporting edges (networks).
	 * 
	 * The additional parameter refRequired determines whether the reference network should always provide support for the consensus edge, or not.
	 * If not, all networks are treated equal. If true, a consensus edge can never be produced if it does not have support in the reference network.
	 * The reference network is determined by trying to cast the input networks to ReferenceNetwork and selecting the one and only unique result.
	 * 
	 * @param networks a set of networks (at least 2).
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param consensusName the name to give to the consensus network
	 * @param ID the ID of the resulting network
	 * @param refRequired whether or not the presence of the edge in the reference network is required for it to be included in the consensus network.
	 * When set to true, this method will raise an IllegalArgumentException when no or more than 1 reference network is found.
	 * @param supportingCutoff the minimal number of networks that need to support an edge in the consensus network
	 * @param weightCutoff the minimal weight that the resulting consensus edges should have to be included
	 * @param minOperator whether or not to take the minimum of the edge weights - if false, the maximum is taken
	 * @param task the task object that keeps track of the progress of this calculation (can be null)
	 * 
	 * @return the consensus network between the input networks. 
	 * 
	 * TODO v3.0: expand this algorithm to be able to deal with n-m node mappings
	 */
	protected ConsensusNetwork calculateConsensusNetwork(Set<Network> networks, TreeEdgeOntology eo, String consensusName, int ID, int supportingCutoff, boolean refRequired, double weightCutoff, boolean minOperator, ScheduledTask task)
	{
		ConsensusNetwork consensus = new ConsensusNetwork(consensusName, ID, networks);

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
			String errormsg = "Please define exactly 1 reference network (" + refNetworks.size() + " found) or change the refRequired parameter to false!";
			throw new IllegalArgumentException(errormsg);
		}

		if (allNetworks.size() == 1)
		{
			Network onlyInput = allNetworks.values().iterator().next();
			consensus.setNodesAndEdges(onlyInput.getNodes(), onlyInput.getEdges());
		}

		Map<String, Set<String>> source2targetIDs = new HashMap<String, Set<String>>();
		Map<String, Set<Node>> nodesByID = new HashMap<String, Set<Node>>();

		for (Network n : networks)
		{
			for (Edge e : n.getEdges())
			{
				Node source = e.getSource();
				String sourceID = source.getID();
				Node target = e.getTarget();
				String targetID = target.getID();

				// add the source-target mapping
				if (!source2targetIDs.containsKey(sourceID))
				{
					source2targetIDs.put(sourceID, new HashSet<String>());
				}
				source2targetIDs.get(sourceID).add(targetID);

				// add the source node by ID
				if (!nodesByID.containsKey(sourceID))
				{
					nodesByID.put(sourceID, new HashSet<Node>());
				}
				nodesByID.get(sourceID).add(source);

				// add the target node by ID
				if (!nodesByID.containsKey(targetID))
				{
					nodesByID.put(targetID, new HashSet<Node>());
				}
				nodesByID.get(targetID).add(target);
			}
		}

		Map<String, Node> allDiffNodes = new HashMap<String, Node>();
		Set<String> roots = eo.retrieveAllSourceRootCats(true);
		
		NodeComparison nc = new NodeComparison();
		EdgeComparison ec = new EdgeComparison(eo);
		EdgeGenerator eg = new EdgeGenerator();
		NetworkCleaning cleaning = new NetworkCleaning(log);

		String networkNames = "";
		for (Network n : networks)
		{
			networkNames += n.getName() + "/";
		}
		networkNames = networkNames.substring(0, networkNames.length() - 1); // cutoff the last "/"
		
		int totalPairs = 0;
		for (String sourceID : source2targetIDs.keySet())
		{
			totalPairs += source2targetIDs.get(sourceID).size();
		}
		
		int ticksPerReport = 1;
		int pairsPerReport = 0;
		
		if (task != null)
		{
			task.setMessage("Calculating consensus network between " + networkNames);
			task.ticksDone(ticksPerReport);
			pairsPerReport = ticksPerReport * totalPairs / task.ticksToGo();
		}

		int pairsProgressed = 0;
		for (String sourceID : source2targetIDs.keySet()) // source node ID 
		{
			Node output_source = nc.getConsensusNode(null, nodesByID.get(sourceID), supportingCutoff);	

			for (String targetID : source2targetIDs.get(sourceID)) // target node ID 
			{
				Node output_target = nc.getConsensusNode(null, nodesByID.get(targetID), supportingCutoff);	
				pairsProgressed++;
				
				// notify the progress listener of our progress
				if (task != null && pairsPerReport > 0 && pairsProgressed % pairsPerReport == 0)
				{
					task.ticksDone(ticksPerReport);
				}
				
				// get all edges from the input network
				Map<String, Map<Integer, Set<EdgeDefinition>>> edgesBySemanticRoot = new HashMap<String, Map<Integer, Set<EdgeDefinition>>>();
				for (String root : roots)
				{
					edgesBySemanticRoot.put(root, new HashMap<Integer, Set<EdgeDefinition>>());
				}
				for (Network n : networks)
				{
					Set<Edge> allEdges = n.getAllEdges(sourceID, targetID);
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

					Map<EdgeDefinition, Set<Integer>> defs = ec.getConsensusEdge(cleanedEdges, supportingNetworks, supportingCutoff, weightCutoff, minOperator);

					// non-void consensus edge
					for (EdgeDefinition def : defs.keySet())
					{
						Set<Integer> supports = defs.get(def);
						if (def.getWeight() > 0)
						{
							if (!allDiffNodes.containsKey(sourceID))
							{
								allDiffNodes.put(sourceID, output_source);
							}
							Node sourceresult = allDiffNodes.get(sourceID);

							if (!allDiffNodes.containsKey(targetID))
							{
								allDiffNodes.put(targetID, output_target);
							}
							Node targetresult = allDiffNodes.get(targetID);

							Set<Condition> conditions = new HashSet<Condition>();
							boolean inReference = false;

							if (supports.isEmpty())
							{
								System.out.println("Found no support for " + sourceresult + " - " + targetresult + " - " + def);
							}

							for (int i : supports)
							{
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
								Edge consensusdiff = new Edge(sourceresult, targetresult, def);
								consensus.addEdge(consensusdiff);
							}
						}
					}
				}
			}
		}
		// notify the progress listener of the fact that we're done (100%)
		if (task != null)
		{
			task.done();
		}

		return consensus;
	}
}
