package be.svlandeg.diffany.usecase.arabidopsis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.StringTokenizer;
import java.util.TreeSet;

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
	public Map<Node, Boolean> getSignificantGenes(OverexpressionData data, double threshold) throws IOException, URISyntaxException
	{
		Map<Node, Boolean> nodes = new HashMap<Node, Boolean>();

		SortedSet<String> ids = data.getArrayIDs();
		SortedSet<String> sign_ids_up = new TreeSet<String>();
		SortedSet<String> sign_ids_down = new TreeSet<String>();
		for (String id : ids)
		{
			double FDR = data.getFDR(id);
			if (FDR <= threshold)
			{
				double FC = data.getFoldchange(id);
				if (FC > 0)
				{
					sign_ids_up.add(id);
				}
				else
				{
					sign_ids_down.add(id);
				}
			}
		}
		for (String id : sign_ids_up)
		{
			nodes.put(new Node(id, id), true);
		}
		for (String id : sign_ids_down)
		{
			nodes.put(new Node(id, id), false);
		}
		return nodes;
	}
	
	/**
	 * TODO
	 * @param targets
	 * @return
	 */
	public Set<Edge> constructVirtualRegulations(Map<Node, Boolean> targets)
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> virtualNodes = new HashMap<String, Node>();
		
		for (Node n : targets.keySet())
		{
			boolean up = targets.get(n);
			String type = "upregulated";
			String regulator = "upregulator";
			if (! up)
			{
				type = "downregulated";
				regulator = "downregulator";
			}
			String fullname = "virtual_" + regulator + "_of_" + n.getID();
			if (! virtualNodes.containsKey(fullname))
			{
				virtualNodes.put(fullname, new Node(fullname, fullname, true));
			}
			Node virtualRegulator = virtualNodes.get(fullname);
			Edge e = new Edge(type, virtualRegulator, n, false);
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
