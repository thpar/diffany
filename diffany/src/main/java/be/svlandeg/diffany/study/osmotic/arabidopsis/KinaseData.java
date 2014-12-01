package be.svlandeg.diffany.study.osmotic.arabidopsis;

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
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.study.osmotic.GenePrinter;

/**
 * This class defines and analyses the data retrieved for Arabidopsis thaliana.
 * 
 * @author Sofie Van Landeghem
 */
public class KinaseData
{
	
	public static final int hubPhos = 30;
	
	private static String kinaseInteractionName = "kinase-targets_20131210.csv";
	private static String kinaseFunctionName = "kinase_activity_go_15102014.tab";
	private static String phosphatName = "phosphat_20130429.csv";
	
	private GenePrinter gp;
	private URI kinaseInteractionLocation;
	private URI kinaseFunctionLocation;
	private URI phosphatLocation;
	
	/**
	 * Constructor: defines the gene printer that can deal with gene synonymy etc.
	 * 
	 * @param gp the gene printer object
	 */
	public KinaseData(GenePrinter gp)
	{
		this.gp = gp;
		
		phosphatLocation = getPhosphat();
		kinaseFunctionLocation = getKinases();
		kinaseInteractionLocation = getKinaseInteractions();
	}

	/**
	 * Retrieve the URI of the Phosphat data
	 * @return the URI of the Phosphat data, or null if the resource could not be located
	 */
	protected URI getPhosphat()
	{
		try
        {
	        return Thread.currentThread().getContextClassLoader().getResource("data/" + phosphatName).toURI();
        }
        catch (URISyntaxException e)
        {
	        System.out.println(" !  Couldn't read " + phosphatName);
        }
		return null;
	}
	
	/**
	 * Retrieve the URI of the kinase activity data
	 * @return the URI of the kinase activity data, or null if the resource could not be located
	 */
	protected URI getKinases()
	{
		try
        {
	        return Thread.currentThread().getContextClassLoader().getResource("data/" + kinaseFunctionName).toURI();
        }
        catch (URISyntaxException e)
        {
	        System.out.println(" !  Couldn't read " + kinaseFunctionName);
        }
		return null;
	}
	
	/**
	 * Retrieve the URI of the kinase interaction data
	 * @return the URI of the kinase interaction data, or null if the resource could not be located
	 */
	protected URI getKinaseInteractions()
	{
		try
        {
	        return Thread.currentThread().getContextClassLoader().getResource("data/" + kinaseInteractionName).toURI();
        }
        catch (URISyntaxException e)
        {
	        System.out.println(" !  Couldn't read " + kinaseInteractionName);
        }
		return null;
	}
	
	/**
	 * Construct a set of kinase-target edges, reading input from a specified URI. This method imposes asymmetry of the read edges.
	 * 
	 * @param includeSelfInteractions whether or not to include self interactions
	 * 
	 * @return the set of PPI edges read from the input file
	 * @throws URISyntaxException when the regulation datafile can not be read properly
	 * @throws IOException when the regulation datafile can not be read properly
	 */
	public Set<Edge> readAllKinaseInteractions(boolean includeSelfInteractions) throws URISyntaxException, IOException
	{
		return readKinaseInteractionsByLocustags(null, null, null, null, includeSelfInteractions);
	}

	/**
	 * Construct a set of kinase interaction edges from a certain input set of nodes, and reading input from a specified URI.
	 * This method can either only include regulations between the nodes themselves, or also include neighbours,
	 * or put a cutoff on minimal number of neighbours to avoid including outliers in the networks.
	 * 
	 * @param source_incl the source nodes (can be null, in which case any node will qualify, except those in source_excl)
	 * @param source_excl the source nodes that are excluded (can be null, in which case any node in source_incl will qualify)
	 * @param target_incl the target nodes (can be null, in which case any node will qualify, except those in target_excl)
	 * @param target_excl the target nodes that are excluded (can be null, in which case any node in target_incl will qualify)
	 * @param includeSelfInteractions whether or not to include self interactions
	 * 
	 * @return the corresponding set of regulation edges
	 * @throws URISyntaxException when the regulation datafile can not be read properly
	 * @throws IOException when the regulation datafile can not be read properly
	 */
	public Set<Edge> readKinaseInteractionsByLocustags(Set<Node> source_incl, Set<Node> source_excl, Set<Node> target_incl, Set<Node> target_excl, boolean includeSelfInteractions) throws URISyntaxException, IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		Map<String, Node> mappedNodes = NodeMapper.getNodesByID(source_incl);
		mappedNodes.putAll(NodeMapper.getNodesByID(source_excl));
		mappedNodes.putAll(NodeMapper.getNodesByID(target_incl));
		mappedNodes.putAll(NodeMapper.getNodesByID(target_excl));

		Set<String> origSourceInclLoci = NodeMapper.getNodeIDs(source_incl);
		Set<String> origSourceExclLoci = NodeMapper.getNodeIDs(source_excl);
		Set<String> origTargetInclLoci = NodeMapper.getNodeIDs(target_incl);
		Set<String> origTargetExclLoci = NodeMapper.getNodeIDs(target_excl);

		BufferedReader reader = new BufferedReader(new FileReader(new File(kinaseInteractionLocation)));

		boolean symmetrical = false;
		Set<String> interactionsRead = new HashSet<String>();

		Map<String, String> mappedTypes = new HashMap<String, String>();
		mappedTypes.put("motif phosphorylation", "phosphorylation");
		mappedTypes.put("peptidearray phosphorylation", "phosphorylation");

		Set<String> excludedTypes = new HashSet<String>();
		excludedTypes.add("interaction");
		excludedTypes.add("pathway");
		excludedTypes.add("regulation");

		String line = reader.readLine();
		line = reader.readLine(); // skip header
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

			boolean include = !(excludedTypes.contains(interaction_type));

			if (include)
			{
				String interactionRead = source_locus + target_locus + interaction_type;

				/* avoid reading the same regulation twice */
				if (!interactionsRead.contains(interactionRead))
				{
					interactionsRead.add(interactionRead);

					boolean foundSource = (source_incl == null || origSourceInclLoci.contains(source_locus)) && (source_excl == null || !origSourceExclLoci.contains(source_locus));
					boolean foundTarget = (target_incl == null || origTargetInclLoci.contains(target_locus)) && (target_excl == null || !origTargetExclLoci.contains(target_locus));

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
	 * @param includePredicted whether or not to include predicted sites
	 * 
	 * @return the corresponding set of locus tags with phosphorylation sites
	 * @throws URISyntaxException when the phosphorylation datafile can not be read properly
	 * @throws IOException when the phosphorylation datafile can not be read properly
	 */
	public Set<String> readPhosphorylationLocusTags(boolean includePredicted) throws URISyntaxException, IOException
	{
		Set<String> locustags = new HashSet<String>();

		BufferedReader reader = new BufferedReader(new FileReader(new File(phosphatLocation)));

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
	 * @return the corresponding set of locus tags with kinase activity
	 * @throws URISyntaxException when the kinase activity datafile can not be read properly
	 * @throws IOException when the kinase activity datafile can not be read properly
	 */
	public Set<String> readKinaseLocusTags() throws URISyntaxException, IOException
	{
		Set<String> locustags = new HashSet<String>();

		BufferedReader reader = new BufferedReader(new FileReader(new File(kinaseFunctionLocation)));
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
