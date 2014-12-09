package be.svlandeg.diffany.study.osmotic.arabidopsis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;

import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.semantics.NodeMapper;

/**
 * This class defines and analyses the regulatory data retrieved for Arabidopsis thaliana.
 * 
 * @author Sofie Van Landeghem
 */
public class RegData
{
	
	private static String atRegName = "reg_net_20100205.tab";  
	private static URI atRegLocation;   
	
	private GenePrinter gp;


	/**
	 * Constructor: defines the gene printer that can deal with gene synonymy etc.
	 * 
	 * @param gp the gene printer object
	 */
	public RegData(GenePrinter gp)
	{
		this.gp = gp;
		atRegLocation = getAtReg();
	}
	
	/**
	 * Retrieve the URI of the regulatory data
	 * @return the URI of the regulatory data, or null if the resource could not be located
	 */
	public URI getAtReg()
	{
		try
        {
	        return Thread.currentThread().getContextClassLoader().getResource("data/arabidopsis/" + atRegName).toURI();
        }
        catch (URISyntaxException e)
        {
	        System.out.println(" !  Couldn't read " + atRegName);
        }
		return null;
	}
	
	/**
	 * Construct a set of regulation edges from a certain input set of nodes, and reading input from a specified URI.
	 * This method can either only include regulations between the nodes themselves, or also include neighbours, or put a cutoff on minimal number of neighbours
	 * to avoid including outliers in the networks.
	 * This method can further remove unspecified regulations, which are not known to be repressors or activators.
	 * 
	 * @param source_nodes the set of input source nodes (or null when any are allowed)
	 * @param target_nodes the set of input target nodes (or null when any are allowed)
	 * @param includeSelfInteractions whether or not to include self interactions
	 * @param includeUnknownPolarities whether or not to include regulatory associations for which we can not determine the polarity (up/down regulation)
	 * 
	 * @return the corresponding set of regulation edges
	 * @throws URISyntaxException when the regulation datafile can not be read properly
	 * @throws IOException when the regulation datafile can not be read properly
	 */
	public Set<Edge> readRegsByLocustags(Set<Node> source_nodes, Set<Node> target_nodes, boolean includeSelfInteractions, boolean includeUnknownPolarities) throws URISyntaxException, IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> mappedNodes = NodeMapper.getNodesByID(source_nodes);
		mappedNodes.putAll(NodeMapper.getNodesByID(target_nodes));

		Set<String> origSourceLoci = NodeMapper.getNodeIDs(source_nodes);
		Set<String> origTargetLoci = NodeMapper.getNodeIDs(target_nodes);

		BufferedReader reader = new BufferedReader(new FileReader(new File(atRegLocation)));

		boolean symmetrical = false;
		Set<String> regsRead = new HashSet<String>();

		String line = reader.readLine();
		line = reader.readLine(); // skip header
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
}
