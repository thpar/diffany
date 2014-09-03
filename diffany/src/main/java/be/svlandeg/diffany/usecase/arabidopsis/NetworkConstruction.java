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
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.semantics.DefaultNodeMapper;
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
	 * @param gp the gene printer object
	 */
	public NetworkConstruction(GenePrinter gp)
	{
		this.gp = gp;
	}
	
	/**
	 * This method will become obsolete once 'expandNetwork' is implemented - delete when not used anymore!
	 * 
	 * @param nodes the overexpressed nodes
	 * @param virtual whether or not to create virtual edges
	 * @param ppi_file the location of the PPI data - or null if you don't want any PPI data
	 * @param reg_file the location of the regulatory data - or null if you don't want any regulatory data
	 * @param selfInteractions whether or not to include self interactions
	 * @param neighbours  whether or not to include direct neighbours
	 * @param includeUnknownReg  whether or not to include unknown regulations
	 * 
	 * @return the set of found edges
	 * @throws IOException when an IO error occurs
	 * @throws URISyntaxException when an input file could not be read
	 * 
	 */
	public Set<Edge> createAllEdgesFromDiffData_old(Map<Node, Double> nodes, boolean virtual, URI ppi_file, URI reg_file, boolean selfInteractions, boolean neighbours, boolean includeUnknownReg) throws URISyntaxException, IOException
	{
		NodeMapper nm = new DefaultNodeMapper();	// TODO: define elsewhere?
		Set<String> origNodes = getNodeIDs(nodes.keySet());
		Set<Edge> edges = new HashSet<Edge>();
		
		if (virtual)
		{
			Set<Edge> virtualEdges = constructVirtualRegulations(nodes);
			System.out.println("  Found " + virtualEdges.size() + " virtual regulations between " + origNodes.size() + " DE genes");
			edges.addAll(virtualEdges);
		}
		
		if (ppi_file != null)
		{
			Set<Edge> PPIedges = readPPIsByLocustags(ppi_file, nodes.keySet(), selfInteractions, neighbours);
			
			Network PPInetwork = new InputNetwork("PPI network", 342, nm);
			PPInetwork.setNodesAndEdges(PPIedges);
			Set<String> expandedNodes = getNodeIDs(PPInetwork.getNodes());
			int expanded = expandedNodes.size();
			expandedNodes.retainAll(origNodes);
			int orig = expandedNodes.size();
			System.out.println("  Found " + PPIedges.size() + " PPIs between " + expanded + " genes, of which " + orig + " DE");
			
			edges.addAll(PPIedges);
		}

		if (reg_file != null)
		{
			// currently, this call does not distinguish between adding neighbours which are targets, and neighbours which are regulators.
			// If this would be important, you could also call the same method a few times to fetch exactly those edges you want.
			Set<Edge> regEdges = readRegsByLocustags(reg_file, nodes.keySet(), nodes.keySet(), selfInteractions, neighbours, neighbours, includeUnknownReg);
			System.out.println("  Found " + regEdges.size() + " regulations");
			edges.addAll(regEdges);
		}
		return edges;
	}
	
	/**
	 * This method defines all the nodes that will be in the Diffany networks to analyse a given set of overexpressed genes.
	 * The set of strict and 'fuzzy' DE genes thus contain all genes that are DE in at least one of the conditions in the experiment.
	 * 
	 * @param nm the object that defines equality between nodes
	 * @param nodes_strict_DE the overexpressed nodes, with stringent criteria (e.g. FDR 0.05)
	 * @param nodes_fuzzy_DE the overexpressed nodes, with less stringent criteria (e.g. FDR 0.1)
	 * @param ppi_file the location of the PPI data
	 * @param reg_file the location of the regulatory data
	 * @param selfInteractions whether or not to include self interactions
	 * @param neighbours  whether or not to include direct neighbours
	 * @param includeUnknownReg  whether or not to include unknown regulations
	 * 
	 * @return the set of found edges
	 * @throws IOException when an IO error occurs
	 * @throws URISyntaxException when an input file could not be read
	 * 
	 */
	public Set<Node> expandNetwork(NodeMapper nm, Set<Node> nodes_strict_DE, Set<Node> nodes_fuzzy_DE, URI ppi_file, URI reg_file, boolean selfInteractions, boolean neighbours, boolean includeUnknownReg) throws URISyntaxException, IOException
	{
		Set<String> origStrictNodesIDs = getNodeIDs(nodes_strict_DE);
		Set<String> origFuzzyNodesIDs = getNodeIDs(nodes_fuzzy_DE);
		
		Set<Node> allNodes = new HashSet<Node>();
		allNodes.addAll(nodes_strict_DE);
		
		Set<Edge> edges = new HashSet<Edge>();
		
		// first expand the node set with PPI neighbours
		if (ppi_file != null)
		{
			Set<Edge> PPIedges = readPPIsByLocustags(ppi_file, nodes_strict_DE, selfInteractions, neighbours);
			
			Network PPInetwork = new InputNetwork("PPI network", 342, nm);
			PPInetwork.setNodesAndEdges(PPIedges);
			Set<String> expandedNodeIDs = getNodeIDs(PPInetwork.getNodes());
			allNodes.addAll(PPInetwork.getNodes());
			
			Set<String> tmp1 = new HashSet<String>(expandedNodeIDs);
			tmp1.retainAll(origStrictNodesIDs);
			int subsetStrict = tmp1.size();
			
			Set<String> tmp2 = new HashSet<String>(expandedNodeIDs);
			tmp2.retainAll(origFuzzyNodesIDs);
			int subsetFuzzy = tmp2.size();
			System.out.println("  Found " + PPIedges.size() + " PPIs between " + expandedNodeIDs.size() + " genes, of which " + subsetStrict + " strict DE and " + subsetFuzzy + " fuzzy DE");
			
			edges.addAll(PPIedges);
		}

		if (reg_file != null)
		{
			// currently, this call does not distinguish between adding neighbours which are targets, and neighbours which are regulators.
			// If this would be important, you could also call the same method a few times to fetch exactly those edges you want.
			Set<Edge> regEdges = readRegsByLocustags(reg_file, nodes_strict_DE, nodes_strict_DE, selfInteractions, neighbours, neighbours, includeUnknownReg);
			System.out.println("  Found " + regEdges.size() + " regulations");
			edges.addAll(regEdges);
		}
		return allNodes;
	}

	/**
	 * Create virtual regulation edges for a collection of targets, which are down-regulated (negative weight) or up-regulated.
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
	 * @param ppi_file the location where to find the tab-delimited PPI data
	 * @param includeSelfInteractions whether or not to include self interactions (e.g. homodimers)
	 * 
	 * @return the set of PPI edges read from the input file
	 * @throws URISyntaxException when the PPI datafile can not be read properly
	 * @throws IOException when the PPI datafile can not be read properly
	 */
	public Set<Edge> readAllPPIs(URI ppi_file, boolean includeSelfInteractions) throws URISyntaxException, IOException
	{
		return readPPIsByLocustags(ppi_file, null, includeSelfInteractions, false, true);
	}

	/**
	 * Construct a set of PPI edges from a certain input set of nodes, and reading input from a specified URI.
	 * This method can either only include PPIs between the nodes themselves, or also include neighbours, or put a cutoff on minimal number of neighbours
	 * to avoid including outliers in the networks. This method imposes symmetry of the read edges.
	 * 
	 * @param ppi_file the location where to find the tab-delimited PPI data
	 * @param nodes the set of input nodes
	 * @param includeSelfInteractions whether or not to include self interactions (e.g. homodimers)
	 * @param includeNeighbours whether or not to include PPI neighbours of the original source set
	 * 
	 * @return the corresponding set of PPI edges
	 * @throws URISyntaxException when the PPI datafile can not be read properly
	 * @throws IOException when the PPI datafile can not be read properly
	 */
	public Set<Edge> readPPIsByLocustags(URI ppi_file, Set<Node> nodes, boolean includeSelfInteractions, boolean includeNeighbours) throws URISyntaxException, IOException
	{
		return readPPIsByLocustags(ppi_file, nodes, includeSelfInteractions, includeNeighbours, false);
	}
	
	/**
	 * Construct a set of PPI edges from a certain input set of nodes, and reading input from a specified URI.
	 * This method can either only include PPIs between the nodes themselves, or also include neighbours, or put a cutoff on minimal number of neighbours
	 * to avoid including outliers in the networks. This method imposes symmetry of the read edges.
	 * 
	 * @param ppi_file the location where to find the tab-delimited PPI data
	 * @param nodes the set of input nodes
	 * @param includeSelfInteractions whether or not to include self interactions (e.g. homodimers)
	 * @param includeNeighbours whether or not to include PPI neighbours of the original source set
	 * 
	 * @return the corresponding set of PPI edges
	 * @throws URISyntaxException when the PPI datafile can not be read properly
	 * @throws IOException when the PPI datafile can not be read properly
	 */
	private Set<Edge> readPPIsByLocustags(URI ppi_file, Set<Node> nodes, boolean includeSelfInteractions, boolean includeNeighbours, boolean includeAll) throws URISyntaxException, IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> mappedNodes = getNodesByID(nodes);
		Set<String> origLoci = getNodeIDs(nodes);

		BufferedReader reader = new BufferedReader(new FileReader(new File(ppi_file)));

		boolean symmetrical = true;
		Set<String> ppisRead = new HashSet<String>();

		String line = reader.readLine();
		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			stok.nextToken();  // String id = 
			String locus1 = stok.nextToken().toLowerCase();
			String locus2 = stok.nextToken().toLowerCase();
			String type = stok.nextToken();

			String ppiRead = locus1 + locus2 + type;
			String ppiReverseRead = locus2 + locus1 + type;

			// avoid reading the same PPI twice
			if (!ppisRead.contains(ppiRead) && !ppisRead.contains(ppiReverseRead))
			{
				ppisRead.add(ppiRead);
				ppisRead.add(ppiReverseRead);

				boolean foundL1 = origLoci.contains(locus1);
				boolean foundL2 = origLoci.contains(locus2);

				// include the interaction when both are in the nodeset, or when one of the two is in the node set and neighbours can be included
				if (includeAll || (foundL1 && foundL2) || (foundL1 && includeNeighbours) || (foundL2 && includeNeighbours))
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
	 * @param reg_file the location where to find the tab-delimited regulatory data
	 * @param source_nodes the set of input source nodes
	 * @param target_nodes the set of input source nodes
	 * @param includeSelfInteractions whether or not to include self interactions
	 * @param includeNeighbourSources whether or not to include neighbour regulators of the original set
	 * @param includeNeighbourTargets whether or not to include neighbour targets of the original set
	 * @param includeUnknownPolarities whether or not to include regulatory associations for which we can not determine the polarity (up/down regulation)
	 * 
	 * @return the corresponding set of regulation edges
	 * @throws URISyntaxException when the regulation datafile can not be read properly
	 * @throws IOException when the regulation datafile can not be read properly
	 */
	public Set<Edge> readRegsByLocustags(URI reg_file, Set<Node> source_nodes, Set<Node> target_nodes, boolean includeSelfInteractions, boolean includeNeighbourSources, boolean includeNeighbourTargets, boolean includeUnknownPolarities) throws URISyntaxException, IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> mappedNodes = getNodesByID(source_nodes);
		mappedNodes.putAll(getNodesByID(target_nodes));
		
		Set<String> origSourceLoci = getNodeIDs(source_nodes);	
		Set<String> origTargetLoci = getNodeIDs(target_nodes);

		BufferedReader reader = new BufferedReader(new FileReader(new File(reg_file)));

		boolean symmetrical = false;
		Set<String> regsRead = new HashSet<String>();

		String line = reader.readLine();
		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			stok.nextToken();	// String source_symbol = 
			String source_locus = stok.nextToken().toLowerCase();
			stok.nextToken();	// String source_family = 
			stok.nextToken();	// String target_symbol = 
			String target_locus = stok.nextToken().toLowerCase();
			stok.nextToken();	// String yesno = 
			stok.nextToken();	// String direct = 
			stok.nextToken();	// String confirmed = 
			stok.nextToken();	// String description = 
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
	
					boolean foundSource = origSourceLoci.contains(source_locus);
					boolean foundTarget = origTargetLoci.contains(target_locus);
	
					// include the interaction when both are in the nodeset
					if ((foundSource && foundTarget) || (foundSource && includeNeighbourTargets) || (foundTarget && includeNeighbourSources))
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

	/**
	 * Retrieve a mapping of the given nodes by their IDs
	 */
	private Map<String, Node> getNodesByID(Set<Node> nodes)
	{
		Map<String, Node> mappedNodes = new HashMap<String, Node>();
		if (nodes != null)
		{
			for (Node n : nodes)
			{
				mappedNodes.put(n.getID(), n);
			}
		}
		return mappedNodes;
	}
	
	/**
	 * Retrieve a unique set of node IDs
	 */
	private Set<String> getNodeIDs(Set<Node> nodes)
	{
		Set<String> IDs = new HashSet<String>();
		if (nodes != null)
		{
			for (Node n : nodes)
			{
				IDs.add(n.getID());
			}
		}
		return IDs;
	}

}
