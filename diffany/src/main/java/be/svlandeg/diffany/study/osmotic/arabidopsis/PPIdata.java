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
 * This class defines and analyses the PPI data retrieved for Arabidopsis thaliana.
 * 
 * @author Sofie Van Landeghem
 */
public class PPIdata
{

	public static final int hubPPI = 10;
	private static String cornetPPIName = "validated_cornet_all_ppi_table_17012012.tab";
	
	private GenePrinter gp;
	private URI cornetPPIlocation;
	
	/**
	 * Constructor: defines the gene printer that can deal with gene synonymy etc.
	 * 
	 * @param gp the gene printer object
	 */
	public PPIdata(GenePrinter gp)
	{
		this.gp = gp;
		cornetPPIlocation = getCornetPPI();
	}

	/**
	 * Retrieve the URI of the CORNET PPI data
	 * @return the URI of the PPI data, or null if the resource could not be located
	 */
	protected URI getCornetPPI()
	{
		try
		{
			return Thread.currentThread().getContextClassLoader().getResource("data/" + cornetPPIName).toURI();
		}
		catch (URISyntaxException e)
		{
			System.out.println(" !  Couldn't read " + cornetPPIName);
		}
		return null;
	}
	
	/**
	 * Construct a set of PPI edges, reading input from a specified URI. This method imposes symmetry of the read edges.
	 * 
	 * @param includeSelfInteractions whether or not to include self interactions (e.g. homodimers)
	 * 
	 * @return the set of PPI edges read from the input file
	 * @throws URISyntaxException when the PPI datafile can not be read properly
	 * @throws IOException when the PPI datafile can not be read properly
	 */
	public Set<Edge> readAllPPIs(boolean includeSelfInteractions) throws URISyntaxException, IOException
	{
		return readPPIsByLocustags(null, null, null, null, includeSelfInteractions);
	}
	

	/**
	 * Construct a set of PPI edges from two input set of nodes, and reading input from a specified URI.
	 * This method can either only include PPIs between the nodes themselves, or also include neighbours (if the second set is null).
	 * This method imposes symmetry of the read edges.
	 * 
	 * @param nodes1_incl the first set of input nodes (can be null, in which case any node will qualify, except those in nodes1_excl)
	 * @param nodes1_excl the nodes that are excluded from the first set of input nodes (can be null, in which case any node in nodes1_incl will qualify)
	 * @param nodes2_incl the second set of input nodes (can be null, in which case any node will qualify, except those in nodes2_excl)
	 * @param nodes2_excl the nodes that are excluded from the second set of input nodes (can be null, in which case any node in nodes2_incl will qualify)
	 * @param includeSelfInteractions whether or not to include self interactions (e.g. homodimers)
	 * 
	 * @return the corresponding set of PPI edges
	 * @throws URISyntaxException when the PPI datafile can not be read properly
	 * @throws IOException when the PPI datafile can not be read properly
	 */
	public Set<Edge> readPPIsByLocustags(Set<Node> nodes1_incl, Set<Node> nodes1_excl, Set<Node> nodes2_incl, Set<Node> nodes2_excl,  boolean includeSelfInteractions) throws URISyntaxException, IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> mappedNodes = NodeMapper.getNodesByID(nodes1_incl);
		mappedNodes.putAll(NodeMapper.getNodesByID(nodes2_incl));
		mappedNodes.putAll(NodeMapper.getNodesByID(nodes1_excl));
		mappedNodes.putAll(NodeMapper.getNodesByID(nodes2_excl));
		Set<String> origLociIncl1 = NodeMapper.getNodeIDs(nodes1_incl);
		Set<String> origLociIncl2 = NodeMapper.getNodeIDs(nodes2_incl);
		Set<String> origLociExcl1 = NodeMapper.getNodeIDs(nodes1_excl);
		Set<String> origLociExcl2 = NodeMapper.getNodeIDs(nodes2_excl);

		BufferedReader reader = new BufferedReader(new FileReader(new File(cornetPPIlocation)));

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

				boolean foundFirstIn1 = (nodes1_incl == null || origLociIncl1.contains(locus1)) && (nodes1_excl == null || ! origLociExcl1.contains(locus1));
				boolean foundFirstIn2 = (nodes2_incl == null || origLociIncl2.contains(locus1)) && (nodes2_excl == null || ! origLociExcl2.contains(locus1));
				boolean foundSecondIn1 = (nodes1_incl == null || origLociIncl1.contains(locus2)) && (nodes1_excl == null || ! origLociExcl1.contains(locus2));
				boolean foundSecondIn2 = (nodes2_incl == null || origLociIncl2.contains(locus2)) && (nodes2_excl == null || ! origLociExcl2.contains(locus2));

				/* include the interaction when both are in one of the nodesets */
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

}
