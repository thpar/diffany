package be.svlandeg.diffany.usecase.arabidopsis;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * This class allows printing gene information in a human-readable format, by taking the original array ID,
 * fetching its corresponding locus tags and Entrez Gene IDs, and gene symbols.
 * 
 * @author Sofie Van Landeghem
 */
public class GenePrinter
{
	
	private Map<String, Set<String>> arrayidmapping;
	private Map<String, String> locusidmapping;
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
		synonymmapping = new MapID().getSynonymMappings(new File(symbolMappingURL.toURI()));
	}
	
	/**
	 * Print an A.th. gene by its array ID
	 * @param arrayID the ID from the A.th microarray dataset
	 */
	public void printGene(String arrayID) 
	{
		Set<String> locusIDs = arrayidmapping.get(arrayID);
		for (String locusID : locusIDs)
		{
			String egid = locusidmapping.get(locusID);
			System.out.print(arrayID + " - " + locusID + " - GID:" + egid);

			System.out.print("\t");
			Set<String> synonyms = new HashSet<String>(synonymmapping.get(egid));
			{
				synonyms.remove(arrayID);
				synonyms.remove(locusID);
				synonyms.remove(egid);
			}
			for (String synonym : synonyms)
			{
				System.out.print(" [" + synonym + "]");
			}
		}
	}

}
