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
	 * @param nodeIDs_strict_DE the overexpressed nodes (not null), with stringent criteria (e.g. FDR 0.05)
	 * @param nodeIDs_fuzzy_DE the overexpressed nodes (may be null), with less stringent criteria (e.g. FDR 0.1, this set does not need to include the strict DE ones)
	 * @param ppi_file the location of the PPI data
	 * @param reg_file the location of the regulatory data (can be null)
	 * @param kinase_file the location of the kinase interaction data (can be null)
	 * @param selfInteractions whether or not to include self interactions
	 * @param neighbours whether or not to include direct neighbours
	 * @param includeUnknownReg whether or not to include unknown regulations
	 * 
	 * @return the set of found edges
	 * @throws IOException when an IO error occurs
	 * @throws URISyntaxException when an input file could not be read
	 */
	public Set<String> expandNetwork(Set<String> nodeIDs_strict_DE, Set<String> nodeIDs_fuzzy_DE, URI ppi_file, URI reg_file, URI kinase_file, boolean selfInteractions, boolean neighbours, boolean includeUnknownReg) throws URISyntaxException, IOException
	{
		// TODO: this method unnecessarily generates Network objects?

		Set<String> allNodeIDs = new HashSet<String>();
		allNodeIDs.addAll(nodeIDs_strict_DE);

		Set<Node> nodes_strict_DE = gp.getNodesByLocusID(nodeIDs_strict_DE);

		/* 1. expand the (strict) DE node set with PPI neighbours */
		if (ppi_file != null)
		{
			Set<Edge> PPIedges = null;
			if (neighbours)
			{
				PPIedges = readPPIsByLocustags(ppi_file, nodes_strict_DE, null, selfInteractions);
			}
			else
			{
				PPIedges = readPPIsByLocustags(ppi_file, nodes_strict_DE, nodes_strict_DE, selfInteractions);
			}

			InputNetwork PPInetwork = new InputNetwork("PPI network", 342, null);
			PPInetwork.setNodesAndEdges(PPIedges);
			Set<String> expandedNodeIDs = NodeMapper.getNodeIDs(PPInetwork.getNodes());
			allNodeIDs.addAll(expandedNodeIDs);
		}

		/* 2. expand the original (strict) DE node set with regulatory and kinase neighbours */
		Set<Edge> regEdges1 = new HashSet<Edge>();
		if (reg_file != null)
		{
			if (neighbours)
			{
				regEdges1.addAll(readRegsByLocustags(reg_file, nodes_strict_DE, null, selfInteractions, includeUnknownReg)); // from our input to their targets
				regEdges1.addAll(readRegsByLocustags(reg_file, null, nodes_strict_DE, selfInteractions, includeUnknownReg)); // from our input to their sources (may result in redundant nodes but these will be cleaned out later)
			}
		}
		if (kinase_file != null)
		{
			if (neighbours)
			{
				regEdges1.addAll(readKinaseInteractionsByLocustags(kinase_file, nodes_strict_DE, null, selfInteractions)); // from our input to their targets
				regEdges1.addAll(readKinaseInteractionsByLocustags(kinase_file, null, nodes_strict_DE, selfInteractions)); // from our input to their sources (may result in redundant nodes but these will be cleaned out later)
			}
		}
		InputNetwork regNetwork1 = new InputNetwork("Regulatory network", 666, null);
		regNetwork1.setNodesAndEdges(regEdges1);
		Set<String> expandedNodeIDs = NodeMapper.getNodeIDs(regNetwork1.getNodes());

		allNodeIDs.addAll(expandedNodeIDs);

		/* 3. add all fuzzy DE nodes which connect to the strict DE nodes or the PPI/regulatory partners */
		if (nodeIDs_fuzzy_DE != null && !nodeIDs_fuzzy_DE.isEmpty())
		{
			Set<Node> allNodes = gp.getNodesByLocusID(allNodeIDs);
			Set<Node> nodes_fuzzy_DE = gp.getNodesByLocusID(nodeIDs_fuzzy_DE);

			Set<Edge> PPIedges = readPPIsByLocustags(ppi_file, allNodes, nodes_fuzzy_DE, selfInteractions);
			InputNetwork PPInetwork = new InputNetwork("PPI network", 342, null);
			PPInetwork.setNodesAndEdges(PPIedges);
			Set<String> expandedPPINodeIDs = NodeMapper.getNodeIDs(PPInetwork.getNodes());

			Set<Edge> regEdges2 = new HashSet<Edge>();
			if (reg_file != null)
			{
				regEdges2.addAll(readRegsByLocustags(reg_file, nodes_fuzzy_DE, allNodes, selfInteractions, includeUnknownReg)); // from fuzzy DE to our combined set
				regEdges2.addAll(readRegsByLocustags(reg_file, allNodes, nodes_fuzzy_DE, selfInteractions, includeUnknownReg)); // from our combined set to fuzzy DE (may result in redundant nodes but these will be cleaned out later)
			}
			if (kinase_file != null)
			{
				regEdges2.addAll(readKinaseInteractionsByLocustags(kinase_file, nodes_fuzzy_DE, allNodes, selfInteractions)); // from fuzzy DE to our combined set
				regEdges2.addAll(readKinaseInteractionsByLocustags(kinase_file, allNodes, nodes_fuzzy_DE, selfInteractions)); // from our combined set to fuzzy DE (may result in redundant nodes but these will be cleaned out later)
			}

			InputNetwork regNetwork2 = new InputNetwork("Regulatory network", 666, null);
			regNetwork2.setNodesAndEdges(regEdges2);
			Set<String> expandedRegNodeIDs = NodeMapper.getNodeIDs(regNetwork2.getNodes());

			allNodeIDs.addAll(expandedPPINodeIDs);
			allNodeIDs.addAll(expandedRegNodeIDs);
		}

		return allNodeIDs;
	}

	/**
	 * Take a copy of the given set of edges and produce new ones that are weighted according to the fold changes of the DE genes.
	 * 
	 * @param eo the edge ontology which can tell which edge types are anti-correlated (e.g. inhibition)
	 * @param nodes the set of nodes that should be used to construct the new edges
	 * @param origEdges the original set of edges (will not be changed!)
	 * @param all_de_nodes all overexpressed nodes (both stringent and less-stringent criteria)
	 * @return a new set of edges, whose weights are changed according to the fold changes of the DE genes
	 */
	public Set<Edge> adjustEdgesByFoldChanges(EdgeOntology eo, Set<Node> nodes, Set<Edge> origEdges, Map<String, Double> all_de_nodes)
	{
		Map<String, Node> mappedNodes = NodeMapper.getNodesByID(nodes);
		
		Set<Edge> resultEdges = new HashSet<Edge>();
		for (Edge e : origEdges)
		{
			String sourceID = e.getSource().getID();
			String targetID = e.getTarget().getID();
			
			double sourceFC = 0;
			if (all_de_nodes.containsKey(sourceID))
			{
				sourceFC = all_de_nodes.get(sourceID);
			}

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
					// this means they are anti-correlated, e.g. by negative regulation
					correlation = false;
				}
			}

			double weight = 1;

			// both are non-DE (we need to keep this if-statement for the following else statements to work properly)
			if (targetFC == 0 && sourceFC == 0)
			{
				weight = 1;
			}

			// both are up-regulated or non-DE: the positive edge weight is (1 + their FC average), or 0 in case they are anti-correlated
			else if (targetFC >= 0 && sourceFC >= 0)
			{
				weight = 1 + ((sourceFC + targetFC) / 2);
				if (!correlation && targetFC != 0) // if target is non-DE, the weight can remain positive
				{
					weight = 0;
				}
			}
			
			// the target is downregulated (below 0), the source is up (above 0) or non-DE -> break the connection unless they are anti-correlated
			else if (targetFC < 0 && sourceFC >= 0)
			{
				weight = 0;
				if (!correlation)
				{
					weight = 1 + ((sourceFC + Math.abs(targetFC)) / 2);
				}
			}

			// the source is downregulated, the target is upregulated or non-DE -> break the connection, even if they are anti-correlated
			// if they were anti-correlated (e.g. inhibition), this is not the case anymore and the differential will be an increase_regulation
			else if (sourceFC < 0 && targetFC >= 0)
			{
				weight = 0;
			}

			// both are down-regulated: the positive edge weight ]0,1[ is their FC average through an exponential function for normalization, 
			// or 0 if they are anti-correlated
			else if (targetFC < 0 && sourceFC < 0)
			{
				weight = 0;
				if (correlation)
				{
					double neg_avg = (sourceFC + targetFC) / 2;

					// this 0.5 factor determines the steepness of the normalization curve
					// after this normalization, the weight will be somewhere between 0 and 1
					weight = Math.exp(0.5 * neg_avg);
				}
			}

			EdgeDefinition newEdge = new EdgeDefinition(e.getDefinition());
			newEdge.setWeight(weight);
			resultEdges.add(new Edge(mappedNodes.get(sourceID), mappedNodes.get(targetID), newEdge));
		}

		return resultEdges;
	}

	/**
	 * Remove the PPI hub biases from a given set of edges. 
	 * To this end, we look for nodes that are only differently connected to one other interaction partner, which is a PPI hub.
	 * In that case, we revert the PPI connection to reference levels (1) because we assume it's not very meaningful, unless both nodes are DE.
	 * 
	 * @param hubs the PPI hubs denoted by their unique IDs
	 * @param nodes the set of nodes that should be used to construct the new edges
	 * @param origEdges the original set of edges (will not be changed!) - should not contain redundancy
	 * @param type the type of connection we want to filter for
	 * @param all_de_nodes the DE genes
	 * @return a new set of edges, which will be a subset of the original set of edges
	 */
	public Set<Edge> filterForHubs(Set<String> hubs, Set<Node> nodes, Set<Edge> origEdges, String type, Set<String> all_de_nodes)
	{
		Map<String, Node> mappedNodes = NodeMapper.getNodesByID(nodes);
		Set<Edge> resultEdges = new HashSet<Edge>();

		Map<String, Integer> diffEdgeCountByNode = new HashMap<String, Integer>();

		/* Count, per node, the number of differentially expressed edges (weights not equal to 1) - of all types */
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
			String sourceID = e.getSource().getID();
			String targetID = e.getTarget().getID();

			int sourceEdgeDiffCount = diffEdgeCountByNode.get(sourceID);
			int targetEdgeDiffCount = diffEdgeCountByNode.get(targetID);

			boolean keepEdge = true;

			/* we only attempt to remove the differential edge when the type is the one we want to filter */
			if (e.getType().equals(type))
			{
				/* the source is only connected by this one differential edge -> candidate for removal (removal means reverting the weight to 1) */
				if (sourceEdgeDiffCount == 1)
				{
					// if target is a hub and the source is not DE
					if (hubs.contains(targetID) && !(all_de_nodes.contains(sourceID)))
					{
						keepEdge = false;
					}
				}
				/* the target is only connected by this one differential edge -> candidate for removal (removal means reverting the weight to 1) */
				if (targetEdgeDiffCount == 1)
				{
					// if source is a hub and the target is not DE
					if (hubs.contains(sourceID) && !(all_de_nodes.contains(targetID)))
					{
						keepEdge = false;
					}
				}
			}
			/* add a copy of the edge to the result set if we want to keep it */
			if (keepEdge)
			{
				Edge copyE = new Edge(mappedNodes.get(sourceID), mappedNodes.get(targetID), new EdgeDefinition(e.getDefinition()));
				resultEdges.add(copyE);
			}
			/* otherwise, put the edge weight to 1 and keep that copy */
			else
			{
				Edge copyE = new Edge(mappedNodes.get(sourceID), mappedNodes.get(targetID), new EdgeDefinition(e.getDefinition()));
				copyE.getDefinition().setWeight(1.0);
				resultEdges.add(copyE);
			}
		}
		return resultEdges;
	}

	/**
	 * Modify the differentially expressed status of a collection of DE nodes, using over/under expression data.
	 * 
	 * @param DEnodeIDs the map of target node IDs with their corresponding fold change values (negative weights refer to down-regulation)
	 * @param nodes the set of existing nodes in the network - they will acquire additional attributes after this method
	 */
	public void modifyDEState(Map<String, Double> DEnodeIDs, Set<Node> nodes)
	{
		String de_attribute = Node.de_attribute;
		for (Node n : nodes)
		{
			String id = n.getID();
			Double FC = DEnodeIDs.get(id);
			if (FC == null || FC == 0)
			{
				n.setAttribute(de_attribute, Node.not_de);
			}
			else if (FC < 0)
			{
				n.setAttribute(de_attribute, Node.downregulated);
			}
			else if (FC > 0)
			{
				n.setAttribute(de_attribute, Node.upregulated);
			}
		}
	}

	/**
	 * Construct a set of PPI edges, reading input from a specified URI. This method imposes symmetry of the read edges.
	 * 
	 * @param ppi_file the location where to find the tab-delimited PPI data
	 * @param includeSelfInteractions whether or not to include self interactions (e.g. homodimers)
	 * 
	 * @return the set of PPI edges read from the input file
	 * @throws URISyntaxException when the PPI datafile can not be read properly
	 * @throws IOException when the PPI datafile can not be read properly
	 */
	public Set<Edge> readAllPPIs(URI ppi_file, boolean includeSelfInteractions) throws URISyntaxException, IOException
	{
		return readPPIsByLocustags(ppi_file, null, null, includeSelfInteractions);
	}

	/**
	 * Construct a set of PPI edges from two input set of nodes, and reading input from a specified URI.
	 * This method can either only include PPIs between the nodes themselves, or also include neighbours (if the second set is null).
	 * This method imposes symmetry of the read edges.
	 * 
	 * @param ppi_file the location where to find the tab-delimited PPI data
	 * @param nodes1 the first set of input nodes (can be null, in which case any node will qualify)
	 * @param nodes2 the second set of input nodes (can be null, in which case any node will qualify)
	 * @param includeSelfInteractions whether or not to include self interactions (e.g. homodimers)
	 * 
	 * @return the corresponding set of PPI edges
	 * @throws URISyntaxException when the PPI datafile can not be read properly
	 * @throws IOException when the PPI datafile can not be read properly
	 */
	public Set<Edge> readPPIsByLocustags(URI ppi_file, Set<Node> nodes1, Set<Node> nodes2, boolean includeSelfInteractions) throws URISyntaxException, IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> mappedNodes = NodeMapper.getNodesByID(nodes1);
		mappedNodes.putAll(NodeMapper.getNodesByID(nodes2));
		Set<String> origLoci1 = NodeMapper.getNodeIDs(nodes1);
		Set<String> origLoci2 = NodeMapper.getNodeIDs(nodes2);

		BufferedReader reader = new BufferedReader(new FileReader(new File(ppi_file)));

		boolean symmetrical = true;
		Set<String> ppisRead = new HashSet<String>();

		String line = reader.readLine();
		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			stok.nextToken();
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

				/* include the interaction when both are in one of the nodesets (this is automatically true for a nodeset which is null) */
				if ((foundFirstIn1 && foundSecondIn2) || (foundSecondIn1 && foundFirstIn2))
				{
					/* include when the loci are different, or when self interactions are allowed */
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
							source = new Node(locus1, symbol);
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
							target = new Node(locus2, symbol);
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
	public Set<Edge> readRegsByLocustags(URI reg_file, Set<Node> source_nodes, Set<Node> target_nodes, boolean includeSelfInteractions, boolean includeUnknownPolarities) throws URISyntaxException, IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> mappedNodes = NodeMapper.getNodesByID(source_nodes);
		mappedNodes.putAll(NodeMapper.getNodesByID(target_nodes));

		Set<String> origSourceLoci = NodeMapper.getNodeIDs(source_nodes);
		Set<String> origTargetLoci = NodeMapper.getNodeIDs(target_nodes);

		BufferedReader reader = new BufferedReader(new FileReader(new File(reg_file)));

		boolean symmetrical = false;
		Set<String> regsRead = new HashSet<String>();

		String line = reader.readLine();
		line = reader.readLine();	// skip header
		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			stok.nextToken(); // source_symbol 
			String source_locus = stok.nextToken().toLowerCase();
			stok.nextToken(); // source_family 
			stok.nextToken(); // target_symbol 
			String target_locus = stok.nextToken().toLowerCase();
			stok.nextToken(); // yesno 
			stok.nextToken(); // direct 
			stok.nextToken(); // confirmed 
			stok.nextToken(); // description 
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

				/* avoid reading the same regulation twice */
				if (!regsRead.contains(regRead))
				{
					regsRead.add(regRead);

					boolean foundSource = (source_nodes == null || origSourceLoci.contains(source_locus));
					boolean foundTarget = (target_nodes == null || origTargetLoci.contains(target_locus));

					/* include the interaction when both are in the nodeset */
					if (foundSource && foundTarget)
					{
						/* include when the loci are different, or when self interactions are allowed */
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
								source = new Node(source_locus, symbol);
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
								target = new Node(target_locus, symbol);
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

	/**
	 * Construct a set of kinase interaction edges from a certain input set of nodes, and reading input from a specified URI.
	 * This method can either only include regulations between the nodes themselves, or also include neighbours, 
	 * or put a cutoff on minimal number of neighbours to avoid including outliers in the networks.
	 * 
	 * @param kinase_interactions_file the location where to find the tab-delimited regulatory data
	 * @param source_nodes the set of input source nodes (or null when any are allowed)
	 * @param target_nodes the set of input target nodes (or null when any are allowed)
	 * @param includeSelfInteractions whether or not to include self interactions
	 * 
	 * @return the corresponding set of regulation edges
	 * @throws URISyntaxException when the regulation datafile can not be read properly
	 * @throws IOException when the regulation datafile can not be read properly
	 */
	public Set<Edge> readKinaseInteractionsByLocustags(URI kinase_interactions_file, Set<Node> source_nodes, Set<Node> target_nodes, boolean includeSelfInteractions) throws URISyntaxException, IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> mappedNodes = NodeMapper.getNodesByID(source_nodes);
		mappedNodes.putAll(NodeMapper.getNodesByID(target_nodes));

		Set<String> origSourceLoci = NodeMapper.getNodeIDs(source_nodes);
		Set<String> origTargetLoci = NodeMapper.getNodeIDs(target_nodes);

		BufferedReader reader = new BufferedReader(new FileReader(new File(kinase_interactions_file)));

		boolean symmetrical = false;
		Set<String> interactionsRead = new HashSet<String>();
		
		Map<String, String> mappedTypes = new HashMap<String, String>();
		mappedTypes.put("motif phosphorylation", "phosphorylation");	
		mappedTypes.put("peptidearray phosphorylation", "phosphorylation");
		
		Set<String> excludedTypes = new HashSet<String>();
		excludedTypes.add("interaction");
		excludedTypes.add("pathway");
		excludedTypes.add("regulation");		// TODO: exclude only when excluding generic regulations
		
		String line = reader.readLine();
		line = reader.readLine(); 		// skip header
		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, ",");
			stok.nextToken(); // Nr
			String type = stok.nextToken().toLowerCase();
			stok.nextToken(); // GO ID
			stok.nextToken(); // GO term
			stok.nextToken(); // MI ID
			stok.nextToken(); // kinase family
			stok.nextToken(); // kin_phos
			String source_locus = stok.nextToken().toLowerCase();
			String target_locus = stok.nextToken().toLowerCase();
			
			String interaction_type = mappedTypes.get(type);
			if (interaction_type == null)
			{
				interaction_type = type;
			}

			boolean include = ! (excludedTypes.contains(interaction_type));

			if (include)
			{
				String interactionRead = source_locus + target_locus + interaction_type;

				/* avoid reading the same regulation twice */
				if (!interactionsRead.contains(interactionRead))
				{
					interactionsRead.add(interactionRead);

					boolean foundSource = (source_nodes == null || origSourceLoci.contains(source_locus));
					boolean foundTarget = (target_nodes == null || origTargetLoci.contains(target_locus));

					/* include the interaction when both are in the nodeset */
					if (foundSource && foundTarget)
					{
						/* include when the loci are different, or when self interactions are allowed */
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
								source = new Node(source_locus, symbol);
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
								target = new Node(target_locus, symbol);
								mappedNodes.put(target_locus, target);
							}
							mappedNodes.put(source_locus, source);
							mappedNodes.put(target_locus, target);
							Edge regulation = new Edge(interaction_type, source, target, symmetrical);
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

	/**
	 * Read a list of (lower-case) locus tags with phosphorylation sites. Either include all, or only those that are experimentally verified.
	 * Those are apparent from the data by the usage of (pS), (pT) or (pY) in the modified peptide string.
	 * 
	 * @param phos_file the location where to find the CSV phosphorylation data
	 * @param includePredicted whether or not to include predicted sites
	 * 
	 * @return the corresponding set of locus tags with phosphorylation sites
	 * @throws URISyntaxException when the phosphorylation datafile can not be read properly
	 * @throws IOException when the phosphorylation datafile can not be read properly
	 */
	public Set<String> readPhosphorylationLocusTags(URI phos_file, boolean includePredicted) throws URISyntaxException, IOException
	{
		Set<String> locustags = new HashSet<String>();

		BufferedReader reader = new BufferedReader(new FileReader(new File(phos_file)));

		String line = reader.readLine();

		// skip header
		line = reader.readLine();

		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, ",");
			String locus = stok.nextToken().toLowerCase();
			stok.nextToken(); // species
			stok.nextToken(); // peptide
			String peptideModified = stok.nextToken();

			// remove the part of the locus tag behind the .
			if (locus.contains("."))
			{
				locus = locus.substring(0, locus.indexOf("."));
			}

			// only try to add this locus tag if we don't have it already
			if (!locustags.contains(locus))
			{
				if (includePredicted || peptideModified.contains("(pS)") || peptideModified.contains("(pT)") || peptideModified.contains("(pY)"))
				{
					locustags.add(locus);
				}
			}
			line = reader.readLine();
		}
		reader.close();

		return locustags;
	}

	/**
	 * Read a list of (lower-case) locus tags with kinase activity. 
	 * The file, downloaded from Gene Ontology, contains all types of annotations and evidence codes, and is currently not filtered further.
	 * 
	 * @param kinase_file the location where to find the tab-delimited kinase activity file
	 * 
	 * @return the corresponding set of locus tags with kinase activity
	 * @throws URISyntaxException when the kinase activity datafile can not be read properly
	 * @throws IOException when the kinase activity datafile can not be read properly
	 */
	public Set<String> readKinaseLocusTags(URI kinase_file) throws URISyntaxException, IOException
	{
		Set<String> locustags = new HashSet<String>();

		BufferedReader reader = new BufferedReader(new FileReader(new File(kinase_file)));
		String line = reader.readLine();

		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			String locus = stok.nextToken().toLowerCase();

			/* Currently, we are not filtering for evidence code as we noticed that quite some known kinases for instance only have the ISS code */
			if (locus.startsWith("at"))
			{
				locustags.add(locus);
			}

			line = reader.readLine();
		}
		reader.close();

		return locustags;
	}

}
