package be.svlandeg.diffany.usecase.arabidopsis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.StringTokenizer;

import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.Node;

/**
 * This class allows to construct networks out of overexpression/coexpression values.
 * 
 * @author Sofie Van Landeghem
 */
public class NetworkConstruction
{

	// TODO make this class more generic / generalizable

	private static String cornetPPIDataFile = "validated_cornet_all_ppi_table_17012012.tab";
	private static String cornetRegDataFile = "reg_net_20100205.tab";

	private GenePrinter gp;

	public NetworkConstruction()
	{
		try
		{
			gp = new GenePrinter();
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	/**
	 * Retrieve all the significant genes in an overexpression dataset, by using the threshold as a minimal cutoff of the FDR values.
	 * @param data the input datasets
	 * @param threshold the FDR cutoff
	 * @return all nodes above the threshold, mapped to their corresponding fold change
	 */
	public Map<Node, Double> getSignificantGenes(OverexpressionData data, double threshold)
	{
		boolean arrayID = data.indexedByRawArrayIDs();

		Map<Node, Double> nodes = new HashMap<Node, Double>();

		SortedSet<String> ids = data.getArrayIDs();
		for (String id : ids)
		{
			double FDR = data.getFDR(id);
			if (FDR <= threshold)
			{
				String symbol = gp.getSymbolByLocusID(id);
				if (arrayID)
				{
					symbol = Arrays.toString(gp.getSymbolByArrayID(id).toArray());
				}
				if (symbol == null)
				{
					symbol = id;
				}
				double FC = data.getFoldchange(id);
				nodes.put(new Node(id.toLowerCase(), symbol, false), FC);
			}
		}
		return nodes;
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
	 * Construct a set of PPI edges from a certain input set of nodes, and reading input from 'cornetPPIDataFile'. 
	 * This method can either only include PPIs between the nodes themselves, or also include neighbours, or put a cutoff on minimal number of neighbours
	 * to avoid including outliers in the networks. This method imposes symmetry of the read edges.
	 * 
	 * @param nodes the set of input nodes
	 * @param includeSelfInteractions whether or not to include self interactions (e.g. homodimers)
	 * @param includeNeighbours whether or not to include PPI neighbours of the original source set
	 * 
	 * @return the corresponding set of PPI edges
	 * @throws URISyntaxException when the PPI datafile from CORNET can not be read properly
	 * @throws IOException when the PPI datafile from CORNET can not be read properly
	 */
	public Set<Edge> readPPIsByLocustags(Set<Node> nodes, boolean includeSelfInteractions, boolean includeNeighbours) throws URISyntaxException, IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> mappedNodes = getNodesByID(nodes);
		Set<String> origLoci = getNodeIDs(nodes);

		URL inputURL = Thread.currentThread().getContextClassLoader().getResource("data/" + cornetPPIDataFile);
		BufferedReader reader = new BufferedReader(new FileReader(new File(inputURL.toURI())));

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
				if ((foundL1 && foundL2) || (foundL1 && includeNeighbours) || (foundL2 && includeNeighbours))
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
	 * Construct a set of regulation edges from a certain input set of nodes, and reading input from 'cornetRegDataFile'. 
	 * This method can either only include regulations between the nodes themselves, or also include neighbours, or put a cutoff on minimal number of neighbours
	 * to avoid including outliers in the networks.
	 * This method can further remove unspecified regulations, which are not known to be repressors or activators.
	 * 
	 * @param source_nodes the set of input source nodes
	 * @param target_nodes the set of input source nodes
	 * @param includeSelfInteractions whether or not to include self interactions
	 * @param includeNeighbourSources whether or not to include neighbour regulators of the original set
	 * @param includeNeighbourTargets whether or not to include neighbour targets of the original set
	 * @param includeUnknownPolarities whether or not to include regulatory associations for which we can not determine the polarity (up/down regulation)
	 * 
	 * @return the corresponding set of regulation edges
	 * @throws URISyntaxException when the regulation datafile from CORNET can not be read properly
	 * @throws IOException when the regulation datafile from CORNET can not be read properly
	 */
	public Set<Edge> readRegsByLocustags(Set<Node> source_nodes, Set<Node> target_nodes, boolean includeSelfInteractions, boolean includeNeighbourSources, boolean includeNeighbourTargets, boolean includeUnknownPolarities) throws URISyntaxException, IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> mappedNodes = getNodesByID(source_nodes);
		mappedNodes.putAll(getNodesByID(target_nodes));
		
		Set<String> origSourceLoci = getNodeIDs(source_nodes);	
		Set<String> origTargetLoci = getNodeIDs(target_nodes);

		URL inputURL = Thread.currentThread().getContextClassLoader().getResource("data/" + cornetRegDataFile);
		BufferedReader reader = new BufferedReader(new FileReader(new File(inputURL.toURI())));

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
		for (Node n : nodes)
		{
			mappedNodes.put(n.getID(), n);
		}
		return mappedNodes;
	}
	
	/**
	 * Retrieve a mapping of the given nodes by their IDs
	 */
	private Set<String> getNodeIDs(Set<Node> nodes)
	{
		Set<String> IDs = new HashSet<String>();
		for (Node n : nodes)
		{
			IDs.add(n.getID());
		}
		return IDs;
	}

}
