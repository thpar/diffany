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
				nodes.put(new Node(id, symbol), FC);
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
			String fullname = regulator + "_of_" + n.getID();
			if (! virtualNodes.containsKey(ID))
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
	public Set<Edge> readPPIsByLocustags(Set<Node> nodes, boolean includeSelfInteractions) throws URISyntaxException, IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> mappedNodes = getNodesByLocustag(nodes);

		URL inputURL = Thread.currentThread().getContextClassLoader().getResource("data/" + cornetDataFile);
		BufferedReader reader = new BufferedReader(new FileReader(new File(inputURL.toURI())));

		String line = reader.readLine();
		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			@SuppressWarnings("unused")
			String id = stok.nextToken();
			String locus1 = stok.nextToken();
			String locus2 = stok.nextToken();
			String type = stok.nextToken();

			if (mappedNodes.keySet().contains(locus1) && mappedNodes.keySet().contains(locus2))
			{
				if (includeSelfInteractions || !locus1.equals(locus2))
				{
					Node source = mappedNodes.get(locus1);
					Node target = mappedNodes.get(locus2);
					boolean symmetrical = true;
					Edge ppi = new Edge(type, source, target, symmetrical);
					edges.add(ppi);
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
	private Map<String, Node> getNodesByLocustag(Set<Node> nodes)
	{
		Map<String, Node> mappedNodes = new HashMap<String, Node>();
		for (Node n : nodes)
		{
			mappedNodes.put(n.getID(), n);
		}
		return mappedNodes;
	}

}
