package be.svlandeg.diffany.usecase.arabidopsis;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.networks.Node;

/**
 * This class allows printing gene information in a human-readable format, by taking the original array ID,
 * fetching its corresponding locus tags and Entrez Gene IDs, and gene symbols.
 * All data is compared in a case independent fashion.
 * 
 * @author Sofie Van Landeghem
 */
public class GenePrinter
{

	private Map<String, Set<String>> arrayidmapping;
	private Map<String, String> locusidmapping;
	private Map<String, String> symbolmapping;
	private Map<String, Set<String>> synonymmapping;


	/**
	 * Create a new GenePrinter, reading the mapping data.
	 * 
	 * @throws IOException when a file can't be read properly
	 * @throws URISyntaxException when a file location can't be parsed properly
	 */
	public GenePrinter() throws IOException, URISyntaxException
	{
		ini();
	}

	/**
	 * Initialize the mapping files from the internal directory structure.
	 * 
	 * @throws IOException when a file can't be read properly
	 * @throws URISyntaxException when a file location can't be parsed properly
	 */
	private void ini() throws IOException, URISyntaxException
	{
		URL arrayMappingURL = Thread.currentThread().getContextClassLoader().getResource("data/affy_ATH1_ID_mapping.tab");
		System.out.println(" Fetching array ID mapping data: " + arrayMappingURL);
		arrayidmapping = new MapID().getAllArrayMappings(new File(arrayMappingURL.toURI()));

		URL locusMappingURL = Thread.currentThread().getContextClassLoader().getResource("data/TAIR10_NCBI_GENEID_mapping.tab");
		System.out.println(" Fetching Locus ID mapping data: " + locusMappingURL);
		locusidmapping = new MapID().getLocusGIDMappings(new File(locusMappingURL.toURI()));

		URL symbolMappingURL = Thread.currentThread().getContextClassLoader().getResource("data/EVEX_synonyms_3702.tab");
		System.out.println(" Fetching gene symbol mapping data: " + symbolMappingURL);
		symbolmapping = new MapID().getSymbolMappings(new File(symbolMappingURL.toURI()));
		synonymmapping = new MapID().getSynonymMappings(new File(symbolMappingURL.toURI()));
	}
	

	/**
	 * Get the synonyms of an A.th. gene by its array ID. There may be multiple entries when an arrayID maps to multiple locus IDs.
	 * 
	 * @param arrayID the ID from the A.th microarray dataset
	 * @return a list of synonyms for this array ID, one concatenated string as entry per locus ID
	 */
	public List<String> getSynonymsByArrayID(String arrayID)
	{
		arrayID = arrayID.toLowerCase();
		List<String> results = new ArrayList<String>();

		Set<String> locusIDs = arrayidmapping.get(arrayID);
		if (locusIDs != null)
		{
			for (String locusID : locusIDs)
			{
				String result = locusID + getSynonymsByLocusID(locusID);
				results.add(result);
			}
		}
		else
		{
			results.add("no_match");
		}
		return results;
	}
	
	/**
	 * Get the official symbol of an A.th. gene by its array ID. There may be multiple entries when an arrayID maps to multiple locus IDs.
	 * 
	 * @param arrayID the ID from the A.th microarray dataset
	 * @return the official symbol, or null if it couldn't be found or this locus ID could not be mapped to Entrez Gene
	 */
	public Set<String> getSymbolByArrayID(String arrayID)
	{
		arrayID = arrayID.toLowerCase();
		Set<String> results = new HashSet<String>();

		Set<String> locusIDs = arrayidmapping.get(arrayID);
		if (locusIDs != null)
		{
			for (String locusID : locusIDs)
			{
				results.add(getSymbolByLocusID(locusID));
			}
		}
		else
		{
			results.add("no_match");
		}
		return results;
	}

	/**
	 * Get the synonyms of an A.th. gene by its locus ID.
	 * 
	 * @param locusID the locus ID of the gene
	 * @return a list of synonyms, concatenated into one string
	 */
	public String getSynonymsByLocusID(String locusID)
	{
		locusID = locusID.toLowerCase();
		String egid = locusidmapping.get(locusID);
		String result = " - GID:" + egid;
		if (egid != null)
		{
			Set<String> synonyms = new HashSet<String>(synonymmapping.get(egid));
			{
				synonyms.remove(locusID);
				synonyms.remove(egid);
			}
			for (String synonym : synonyms)
			{
				result += " [" + synonym + "]";
			}
		}
		return result;
	}
	
	/**
	 * Get the official symbol of an A.th. gene by its locus ID.
	 * 
	 * @param locusID the locus ID of the gene
	 * @return the official symbol, or null if it couldn't be found or this locus ID could not be mapped to Entrez Gene
	 */
	public String getSymbolByLocusID(String locusID)
	{
		locusID = locusID.toLowerCase();
		String egid = locusidmapping.get(locusID);
		if (egid != null)
		{
			return symbolmapping.get(egid);
		}
		return null;
	}
	
	/**
	 * Retrieve a set of nodes from their node IDs, by using the GenePrinter to search for their symbol
	 * @param nodeIDs the locus IDs of the nodes we want to generate
	 * @return the set of nodes, one for each unique locus ID, with their symbol determined by this gene printer
	 */
	public Set<Node> getNodesByLocusID(Set<String> nodeIDs)
	{
		Set<Node> nodes = new HashSet<Node>();
		for (String locusID : nodeIDs)
		{
			String symbol = getSymbolByLocusID(locusID);
			if (symbol == null)
			{
				symbol = locusID;
			}
			nodes.add(new Node(locusID, symbol, false));
		}
		return nodes;
	}

}
