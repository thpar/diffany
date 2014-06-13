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

	private static String cornetDataFile = "validated_cornet_all_ppi_table_17012012.tab";

	/**
	 * 
	 * @param datasets
	 * @throws URISyntaxException 
	 * @throws IOException 
	 */
	public Map<Node, Double> getSignificantGenes(OverexpressionData data, double threshold) throws IOException, URISyntaxException
	{
		GenePrinter gp = new GenePrinter();

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
				nodes.put(new Node(id.toLowerCase(), symbol), FC);
			}
		}
		return nodes;
	}

	/**
	 * TODO
	 * @param targets
	 * @return
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
	 * TODO
	 * @param locusIDs
	 * @return
	 * @throws URISyntaxException
	 * @throws IOException
	 */
	public Set<Edge> readPPIsByLocustags(Set<Node> nodes, boolean includeSelfInteractions, boolean includeNeighbours, int min_neighbourcount) throws URISyntaxException, IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> origNodes = getNodesByID(nodes);
		Map<String, Node> mappedNodes = new HashMap<String, Node>();

		URL inputURL = Thread.currentThread().getContextClassLoader().getResource("data/" + cornetDataFile);
		BufferedReader reader = new BufferedReader(new FileReader(new File(inputURL.toURI())));

		boolean symmetrical = true;
		Set<String> ppisRead = new HashSet<String>();
		
		// for each new neighbour, store its interaction partners from the original set of nodes, to check whether there are more than min_neighbourcount
		Map<String, Set<String>> neighbourCounts = new HashMap<String, Set<String>>();

		String line = reader.readLine();
		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			@SuppressWarnings("unused")
			String id = stok.nextToken();
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
				
				boolean foundL1 = origNodes.keySet().contains(locus1);
				boolean foundL2 = origNodes.keySet().contains(locus2);

				// include the interaction when both are in the nodeset
				if (foundL1 && foundL2)
				{
					// include when the loci are different, or when self interactions are allowed
					if (includeSelfInteractions || !locus1.equals(locus2))
					{
						Node source = origNodes.get(locus1);
						Node target = origNodes.get(locus2);
						mappedNodes.put(locus1, source);
						mappedNodes.put(locus2, target);
						Edge ppi = new Edge(type, source, target, symmetrical);
						edges.add(ppi);
					}
				}
				if (foundL1 && !foundL2)
				{
					if (! neighbourCounts.containsKey(locus2))
					{
						neighbourCounts.put(locus2, new HashSet<String>());
					}
					neighbourCounts.get(locus2).add(locus1);
				}
				
				if (foundL2 && !foundL1)
				{
					if (! neighbourCounts.containsKey(locus1))
					{
						neighbourCounts.put(locus1, new HashSet<String>());
					}
					neighbourCounts.get(locus1).add(locus2);
				}
			}
			line = reader.readLine();
		}
		reader.close();
		
		// Check which neighbours can be added to the network (when they are sufficiently connected)
		Set<String> allowedNeighbours = new HashSet<String>();
		for (String neighbour : neighbourCounts.keySet())
		{
			if (neighbourCounts.get(neighbour).size() >= min_neighbourcount)
			{
				allowedNeighbours.add(neighbour);
			}
		}
		
		// read a second time to include the proper neighbours
		reader = new BufferedReader(new FileReader(new File(inputURL.toURI())));

		ppisRead = new HashSet<String>();

		line = reader.readLine();
		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			@SuppressWarnings("unused")
			String id = stok.nextToken();
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
				
				boolean foundL1 = origNodes.keySet().contains(locus1);
				boolean foundL2 = origNodes.keySet().contains(locus2);
				
				boolean neighbourL1 = allowedNeighbours.contains(locus1);
				boolean neighbourL2 = allowedNeighbours.contains(locus2);

				// include the interaction when both are in the nodeset, or when neighbours can be included and either one is in the node set
				if ((foundL1 && !foundL2 && neighbourL2) || (foundL2 && !foundL1 && neighbourL1))
				{
					// include when the loci are different, or when self interactions are allowed
					if (includeSelfInteractions || !locus1.equals(locus2))
					{
						Node source = mappedNodes.get(locus1);
						if (source == null)
						{
							source = new Node(locus1);
							mappedNodes.put(locus1, source);
						}
						
						Node target = mappedNodes.get(locus2);
						if (target == null)
						{
							target = new Node(locus2);
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
	 * TODO
	 * @param nodes
	 * @return
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

}
