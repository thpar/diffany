package be.svlandeg.diffany.usecase.arabidopsis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;

import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.EdgeDefinition;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.semantics.EdgeOntology;
import be.svlandeg.diffany.core.semantics.NodeMapper;

/**
 * This class allows to construct networks out of overexpression/coexpression values.
 * 
 * @author Sofie Van Landeghem
 */
public class NetworkConstruction
{

	// TODO v2.1 make this class more generic / generalizable

	private GenePrinter gp;

	/**
	 * Constructor: defines the gene printer that can deal with gene synonymy etc.
	 * 
	 * @param gp the gene printer object
	 */
	public NetworkConstruction(GenePrinter gp)
	{
		this.gp = gp;
	}

	/**
	 * This method defines all the nodes that will be in the Diffany networks to analyse a given set of overexpressed genes.
	 * The set of strict and 'fuzzy' DE genes thus contain all genes that are DE in at least one of the conditions in the experiment.
	 * 
	 * @param nm the object that defines equality between nodes
	 * @param nodeIDs_strict_DE the overexpressed nodes (not null), with stringent criteria (e.g. FDR 0.05)
	 * @param nodeIDs_fuzzy_DE the overexpressed nodes (may be null), with less stringent criteria (e.g. FDR 0.1, this set does not need to include the strict DE ones)
	 * @param ppi_file the location of the PPI data
	 * @param reg_file the location of the regulatory data
	 * @param selfInteractions whether or not to include self interactions
	 * @param neighbours whether or not to include direct neighbours
	 * @param includeUnknownReg whether or not to include unknown regulations
	 * 
	 * @return the set of found edges
	 * @throws IOException when an IO error occurs
	 * @throws URISyntaxException when an input file could not be read
	 */
	public Set<String> expandNetwork(NodeMapper nm, Set<String> nodeIDs_strict_DE, Set<String> nodeIDs_fuzzy_DE, URI ppi_file, URI reg_file, boolean selfInteractions, boolean neighbours, boolean includeUnknownReg) throws URISyntaxException, IOException
	{
		Set<String> allNodeIDs = new HashSet<String>();
		allNodeIDs.addAll(nodeIDs_strict_DE);

		Set<Node> nodes_strict_DE = gp.getNodesByLocusID(nodeIDs_strict_DE);

		// first expand the (strict) DE node set with PPI neighbours
		if (ppi_file != null)
		{
			Set<Edge> PPIedges = null;
			if (neighbours)
			{
				PPIedges = readPPIsByLocustags(nm, ppi_file, nodes_strict_DE, null, selfInteractions);
			}
			else
			{
				PPIedges = readPPIsByLocustags(nm, ppi_file, nodes_strict_DE, nodes_strict_DE, selfInteractions);
			}

			InputNetwork PPInetwork = new InputNetwork("PPI network", 342, nm);
			PPInetwork.setNodesAndEdges(PPIedges);
			Set<String> expandedNodeIDs = nm.getNodeIDs(PPInetwork.getNodes());
			allNodeIDs.addAll(expandedNodeIDs);
		}

		// then expand the original (strict) DE node set with regulatory neighbours
		if (reg_file != null)
		{

			Set<Edge> regEdges = null;
			if (neighbours)
			{
				regEdges = readRegsByLocustags(nm, reg_file, nodes_strict_DE, null, selfInteractions, includeUnknownReg); // from our input to their targets
				regEdges.addAll(readRegsByLocustags(nm, reg_file, null, nodes_strict_DE, selfInteractions, includeUnknownReg)); // from our input to their sources (may result in redundant nodes but these will be cleaned out later)
			}

			InputNetwork regNetwork = new InputNetwork("Regulatory network", 666, nm);
			regNetwork.setNodesAndEdges(regEdges);
			Set<String> expandedNodeIDs = nm.getNodeIDs(regNetwork.getNodes());
			allNodeIDs.addAll(expandedNodeIDs);
		}

		// finally, add all fuzzy DE nodes which connect to the strict DE nodes or the PPI/regulatory partners
		if (nodeIDs_fuzzy_DE != null && !nodeIDs_fuzzy_DE.isEmpty())
		{
			Set<Node> allNodes = gp.getNodesByLocusID(allNodeIDs);
			Set<Node> nodes_fuzzy_DE = gp.getNodesByLocusID(nodeIDs_fuzzy_DE);
			Set<Edge> PPIedges = readPPIsByLocustags(nm, ppi_file, allNodes, nodes_fuzzy_DE, selfInteractions);
			InputNetwork PPInetwork = new InputNetwork("PPI network", 342, nm);
			PPInetwork.setNodesAndEdges(PPIedges);
			Set<String> expandedPPINodeIDs = nm.getNodeIDs(PPInetwork.getNodes());

			Set<Edge> regEdges = readRegsByLocustags(nm, reg_file, nodes_fuzzy_DE, allNodes, selfInteractions, includeUnknownReg); // from fuzzy DE to our combined set
			regEdges.addAll(readRegsByLocustags(nm, reg_file, allNodes, nodes_fuzzy_DE, selfInteractions, includeUnknownReg)); // from our combined set to fuzzy DE (may result in redundant nodes but these will be cleaned out later)

			InputNetwork regNetwork = new InputNetwork("Regulatory network", 666, nm);
			regNetwork.setNodesAndEdges(regEdges);
			Set<String> expandedRegNodeIDs = nm.getNodeIDs(regNetwork.getNodes());
			allNodeIDs.addAll(expandedRegNodeIDs);

			allNodeIDs.addAll(expandedPPINodeIDs);
			allNodeIDs.addAll(expandedRegNodeIDs);
		}

		return allNodeIDs;
	}

	/**
	 * Take a copy of the given set of edges and produce new ones that are weighted according to the fold changes of the DE genes.
	 * 
	 * @param eo the edge ontology which can tell which edge types are anti-correlated (e.g. inhibition)
	 * @param origEdges the original set of edges (will not be changed!)
	 * @param all_de_nodes all overexpressed nodes (both stringent and less-stringent criteria)
	 * @return a new set of edges, whose weights are changed according to the fold changes of the DE genes
	 */
	public Set<Edge> adjustEdgesByFoldChanges(EdgeOntology eo, Set<Edge> origEdges, Map<String, Double> all_de_nodes)
	{
		Set<Edge> resultEdges = new HashSet<Edge>();
		for (Edge e : origEdges)
		{
			Edge newEdge = new Edge(e.getSource(), e.getTarget(), new EdgeDefinition(e.getDefinition()));

			String sourceID = newEdge.getSource().getID();
			double sourceFC = 0;
			if (all_de_nodes.containsKey(sourceID))
			{
				sourceFC = all_de_nodes.get(sourceID);
			}

			String targetID = newEdge.getTarget().getID();
			double targetFC = 0;
			if (all_de_nodes.containsKey(targetID))
			{
				targetFC = all_de_nodes.get(targetID);
			}

			boolean correlation = true;
			if (!e.isSymmetrical())
			{
				String edgeType = e.getType();
				String edgeCat = eo.getSourceCategory(edgeType);
				if (eo.getAllNegSourceCategories().contains(edgeCat))
				{
					correlation = false; // this means they are anti-correlated, e.g. by negative regulation
				}
			}

			double weight = 1;

			// both are non-DE (we need to keep this if-statement for the following else statements to work properly)
			if (targetFC == 0 && sourceFC == 0)
			{
				weight = 1;
			}

			// both are up-regulated or equal: the positive edge weight is (1 + their FC average), or 0 in case they are anti-correlated
			else if (targetFC >= 0 && sourceFC >= 0)
			{
				weight = 1 + ((sourceFC + targetFC) / 2);
				if (!correlation && targetFC != 0) // if target is equal, the weight can remain positive
				{
					weight = 0;
				}
			}

			// source is downregulated but target is equal: 0 weight
			else if (targetFC == 0 && sourceFC < 0)
			{
				weight = 0;
			}

			// both are down-regulated: the positive edge weight is their FC average through an exponential function for normalization, 
			// or 0 if they are anti-correlated
			else if (targetFC < 0 && sourceFC < 0)
			{
				weight = 0;
				if (correlation)
				{
					double neg_avg = (sourceFC + targetFC) / 2;

					// this 0.5 factor determines the steepness of the normalization curve
					weight = Math.exp(0.5 * neg_avg);
				}
			}

			// one of the two is downregulated (below 0), the other up (above 0) -> break the connection unless they are anti-correlated
			else if ((targetFC < 0 && sourceFC >= 0) || (targetFC > 0 && sourceFC < 0))
			{
				weight = 0;
				if (!correlation)
				{
					weight = 1 + ((Math.abs(sourceFC) + Math.abs(targetFC)) / 2);
				}
			}

			newEdge.getDefinition().setWeight(weight);
			resultEdges.add(newEdge);
		}

		return resultEdges;
	}

	/**
	 * Remove the PPI hub biases from a given set of edges. To this end, we look for nodes that are only differently connected to one other interaction partner, which is a PPI hub.
	 * In that case, we revert the PPI connection to reference levels (1) because we assume it's not very meaningful, unless both nodes are DE.
	 * 
	 * @param hubs the PPI hubs denoted by their unique IDs
	 * @param origEdges the original set of edges (will not be changed!)
	 * @param type the type of connection we want to filter for
	 * @param all_de_nodes the DE genes
	 * @return a new set of edges, which will be a subset of the original set of edges
	 */
	public Set<Edge> filterForHubs(Set<String> hubs, Set<Edge> origEdges, String type, Set<String> all_de_nodes)
	{
		Set<Edge> resultEdges = new HashSet<Edge>();

		Map<String, Integer> diffEdgeCountByNode = new HashMap<String, Integer>();

		// Count, per node, the number of differentially expressed edges (weights not equal to 1) - of all types
		for (Edge e : origEdges)
		{
			String sourceID = e.getSource().getID();
			String targetID = e.getTarget().getID();

			if (!diffEdgeCountByNode.containsKey(sourceID))
			{
				diffEdgeCountByNode.put(sourceID, 0);
			}
			if (!diffEdgeCountByNode.containsKey(targetID))
			{
				diffEdgeCountByNode.put(targetID, 0);
			}

			if (e.getWeight() != 1.0)
			{
				diffEdgeCountByNode.put(sourceID, diffEdgeCountByNode.get(sourceID) + 1);
				diffEdgeCountByNode.put(targetID, diffEdgeCountByNode.get(targetID) + 1);
			}
		}

		for (Edge e : origEdges)
		{
			Node source = e.getSource();
			Node target = e.getTarget();

			String sourceID = source.getID();
			String targetID = target.getID();

			int sourceEdgeDiffCount = diffEdgeCountByNode.get(sourceID);
			int targetEdgeDiffCount = diffEdgeCountByNode.get(targetID);

			boolean keepEdge = true;

			// we only attempt to remove the differential edge when the type if the one we want to filter
			if (e.getType().equals(type))
			{
				// either source or target is only connected by this one differential edge -> candidate for removal (removal means reverting the weight to 1)
				if (sourceEdgeDiffCount == 1 || targetEdgeDiffCount == 1)
				{
					// if either source or target is a hub, we will probably not keep the edge
					if (hubs.contains(sourceID) || hubs.contains(targetID))
					{
						// unless both source and target are DE, then we do keep the edge
						if (!(all_de_nodes.contains(sourceID) && all_de_nodes.contains(targetID)))
						{
							keepEdge = false;
						}
					}
				}
			}
			// add a copy of the edge to the result set if we want to keep it
			if (keepEdge)
			{
				Edge copyE = new Edge(source, target, new EdgeDefinition(e.getDefinition()));
				resultEdges.add(copyE);
			}
			// otherwise, put the edge weight to 1 and keep that copy
			else
			{
				Edge copyE = new Edge(source, target, new EdgeDefinition(e.getDefinition()));
				copyE.getDefinition().setWeight(1.0);
				resultEdges.add(copyE);
			}
		}
		return resultEdges;
	}

	/**
	 * Create virtual regulation edges for a collection of targets, which are down-regulated (negative weight) or up-regulated.
	 * This method is currently not used anymore, as fold changes are encoded into the PPI and regulatory edges.
	 * 
	 * @param targets the map of target nodes with their corresponding regulation weights (negative weights refer to down-regulation)
	 * @return the set of corresponding virtual edges
	 */
	public Set<Edge> constructVirtualRegulations(Map<Node, Double> targets)
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> virtualNodes = new HashMap<String, Node>();

		for (Node n : targets.keySet())
		{
			double FC = targets.get(n);

			String type = "upregulated";
			String regulator = "upregulator";

			if (FC < 0)
			{
				type = "downregulated";
				regulator = "downregulator";
			}

			String ID = type.charAt(0) + "_" + n.getID();
			String fullname = regulator + "_of_" + n.getDisplayName();
			if (!virtualNodes.containsKey(ID))
			{
				virtualNodes.put(ID, new Node(ID, fullname, true));
			}
			Node virtualRegulator = virtualNodes.get(ID);
			Edge e = new Edge(type, virtualRegulator, n, false, Math.abs(FC), false);
			edges.add(e);
		}

		return edges;
	}

	/**
	 * Construct a set of PPI edges, reading input from a specified URI. This method imposes symmetry of the read edges.
	 * 
	 * @param nm the object that defines equality between nodes
	 * @param ppi_file the location where to find the tab-delimited PPI data
	 * @param includeSelfInteractions whether or not to include self interactions (e.g. homodimers)
	 * 
	 * @return the set of PPI edges read from the input file
	 * @throws URISyntaxException when the PPI datafile can not be read properly
	 * @throws IOException when the PPI datafile can not be read properly
	 */
	public Set<Edge> readAllPPIs(NodeMapper nm, URI ppi_file, boolean includeSelfInteractions) throws URISyntaxException, IOException
	{
		return readPPIsByLocustags(nm, ppi_file, null, null, includeSelfInteractions);
	}

	/**
	 * Construct a set of PPI edges from two input set of nodes, and reading input from a specified URI.
	 * This method can either only include PPIs between the nodes themselves, or also include neighbours (if the second set is null).
	 * This method imposes symmetry of the read edges.
	 * 
	 * @param nm the object that defines equality between nodes
	 * @param ppi_file the location where to find the tab-delimited PPI data
	 * @param nodes1 the first set of input nodes (can be null, in which case any node will qualify)
	 * @param nodes2 the second set of input nodes (can be null, in which case any node will qualify)
	 * @param includeSelfInteractions whether or not to include self interactions (e.g. homodimers)
	 * 
	 * @return the corresponding set of PPI edges
	 * @throws URISyntaxException when the PPI datafile can not be read properly
	 * @throws IOException when the PPI datafile can not be read properly
	 */
	public Set<Edge> readPPIsByLocustags(NodeMapper nm, URI ppi_file, Set<Node> nodes1, Set<Node> nodes2, boolean includeSelfInteractions) throws URISyntaxException, IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> mappedNodes = nm.getNodesByID(nodes1);
		mappedNodes.putAll(nm.getNodesByID(nodes2));
		Set<String> origLoci1 = nm.getNodeIDs(nodes1);
		Set<String> origLoci2 = nm.getNodeIDs(nodes2);

		BufferedReader reader = new BufferedReader(new FileReader(new File(ppi_file)));

		boolean symmetrical = true;
		Set<String> ppisRead = new HashSet<String>();

		String line = reader.readLine();
		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			stok.nextToken(); // String id = 
			String locus1 = stok.nextToken().toLowerCase();
			String locus2 = stok.nextToken().toLowerCase();
			String type = stok.nextToken();
			String dataSource = stok.nextToken();

			boolean unconfirmed = dataSource.contains("non-confirmed");

			String ppiRead = locus1 + locus2 + type;
			String ppiReverseRead = locus2 + locus1 + type;

			// avoid reading the same PPI twice
			if (!ppisRead.contains(ppiRead) && !ppisRead.contains(ppiReverseRead) && !unconfirmed)
			{
				ppisRead.add(ppiRead);
				ppisRead.add(ppiReverseRead);

				boolean foundFirstIn1 = (nodes1 == null || origLoci1.contains(locus1));
				boolean foundFirstIn2 = (nodes2 == null || origLoci2.contains(locus1));
				boolean foundSecondIn1 = (nodes1 == null || origLoci1.contains(locus2));
				boolean foundSecondIn2 = (nodes2 == null || origLoci2.contains(locus2));

				// include the interaction when both are in one of the nodesets (this is automatically true for a nodeset which is null)
				if ((foundFirstIn1 && foundSecondIn2) || (foundSecondIn1 && foundFirstIn2))
				{
					// include when the loci are different, or when self interactions are allowed
					if (includeSelfInteractions || !locus1.equals(locus2))
					{
						Node source = mappedNodes.get(locus1);
						if (source == null)
						{
							String symbol = gp.getSymbolByLocusID(locus1);
							if (symbol == null)
							{
								symbol = locus1;
							}
							source = new Node(locus1, symbol, false);
							mappedNodes.put(locus1, source);
						}

						Node target = mappedNodes.get(locus2);
						if (target == null)
						{
							String symbol = gp.getSymbolByLocusID(locus2);
							if (symbol == null)
							{
								symbol = locus2;
							}
							target = new Node(locus2, symbol, false);
							mappedNodes.put(locus2, target);
						}
						Edge ppi = new Edge(type, source, target, symmetrical);
						edges.add(ppi);
					}
				}
			}
			line = reader.readLine();
		}
		reader.close();

		return edges;
	}

	/**
	 * Construct a set of regulation edges from a certain input set of nodes, and reading input from a specified URI.
	 * This method can either only include regulations between the nodes themselves, or also include neighbours, or put a cutoff on minimal number of neighbours
	 * to avoid including outliers in the networks.
	 * This method can further remove unspecified regulations, which are not known to be repressors or activators.
	 * 
	 * @param nm the object that defines equality between nodes
	 * @param reg_file the location where to find the tab-delimited regulatory data
	 * @param source_nodes the set of input source nodes (or null when any are allowed)
	 * @param target_nodes the set of input target nodes (or null when any are allowed)
	 * @param includeSelfInteractions whether or not to include self interactions
	 * @param includeUnknownPolarities whether or not to include regulatory associations for which we can not determine the polarity (up/down regulation)
	 * 
	 * @return the corresponding set of regulation edges
	 * @throws URISyntaxException when the regulation datafile can not be read properly
	 * @throws IOException when the regulation datafile can not be read properly
	 */
	public Set<Edge> readRegsByLocustags(NodeMapper nm, URI reg_file, Set<Node> source_nodes, Set<Node> target_nodes, boolean includeSelfInteractions, boolean includeUnknownPolarities) throws URISyntaxException, IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> mappedNodes = nm.getNodesByID(source_nodes);
		mappedNodes.putAll(nm.getNodesByID(target_nodes));

		Set<String> origSourceLoci = nm.getNodeIDs(source_nodes);
		Set<String> origTargetLoci = nm.getNodeIDs(target_nodes);

		BufferedReader reader = new BufferedReader(new FileReader(new File(reg_file)));

		boolean symmetrical = false;
		Set<String> regsRead = new HashSet<String>();

		String line = reader.readLine();
		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			stok.nextToken(); // String source_symbol = 
			String source_locus = stok.nextToken().toLowerCase();
			stok.nextToken(); // String source_family = 
			stok.nextToken(); // String target_symbol = 
			String target_locus = stok.nextToken().toLowerCase();
			stok.nextToken(); // String yesno = 
			stok.nextToken(); // String direct = 
			stok.nextToken(); // String confirmed = 
			stok.nextToken(); // String description = 
			String type = stok.nextToken().toLowerCase();

			boolean include = true;
			if (type.equals("unknown"))
			{
				type = "unknown_regulation";
				include = includeUnknownPolarities;
			}

			if (include)
			{
				String regRead = source_locus + target_locus + type;

				// avoid reading the same regulation twice
				if (!regsRead.contains(regRead))
				{
					regsRead.add(regRead);

					boolean foundSource = (source_nodes == null || origSourceLoci.contains(source_locus));
					boolean foundTarget = (target_nodes == null || origTargetLoci.contains(target_locus));

					// include the interaction when both are in the nodeset
					if (foundSource && foundTarget)
					{
						// include when the loci are different, or when self interactions are allowed
						if (includeSelfInteractions || !source_locus.equals(target_locus))
						{
							Node source = mappedNodes.get(source_locus);
							if (source == null)
							{
								String symbol = gp.getSymbolByLocusID(source_locus);
								if (symbol == null)
								{
									symbol = source_locus;
								}
								source = new Node(source_locus, symbol, false);
								mappedNodes.put(source_locus, source);
							}
							Node target = mappedNodes.get(target_locus);
							if (target == null)
							{
								String symbol = gp.getSymbolByLocusID(target_locus);
								if (symbol == null)
								{
									symbol = target_locus;
								}
								target = new Node(target_locus, symbol, false);
								mappedNodes.put(target_locus, target);
							}
							mappedNodes.put(source_locus, source);
							mappedNodes.put(target_locus, target);
							Edge regulation = new Edge(type, source, target, symmetrical);
							edges.add(regulation);
						}
					}
				}
			}
			line = reader.readLine();
		}
		reader.close();

		return edges;
	}

}
