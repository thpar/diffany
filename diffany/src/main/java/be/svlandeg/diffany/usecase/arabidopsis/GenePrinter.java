package be.svlandeg.diffany.usecase.arabidopsis;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * This class allows printing gene information in a human-readable format
 * 
 * @author Sofie Van Landeghem
 */
public class GenePrinter
{
	
	private Map<String, Set<String>> arrayidmapping;
	private Map<String, String> locusidmapping;
	private Map<String, String> symbolmapping;
	private Map<String, Set<String>> synonymmapping;
	
	public GenePrinter() throws IOException, URISyntaxException 
	{
		ini();
	}
	
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
	 * Print an A.th. gene by its array ID
	 * @throws URISyntaxException 
	 * @throws IOException 
	 */
	public void printGene(String arrayID) 
	{
		Set<String> locusIDs = arrayidmapping.get(arrayID);
		for (String locusID : locusIDs)
		{
			String egid = locusidmapping.get(locusID);
			String symbol = symbolmapping.get(egid);
			System.out.print(arrayID + " - " + locusID + " - GID:" + egid);
			/*if (symbol != null && !symbol.equals(locusID))
			{
				System.out.print(" - " + symbol);
			}*/
			System.out.print("\t");
			Set<String> synonyms = new HashSet<String>(synonymmapping.get(egid));
			{
				synonyms.remove(arrayID);
				synonyms.remove(locusID);
				synonyms.remove(egid);
				synonyms.remove(symbol);
			}
			for (String synonym : synonyms)
			{
				System.out.print(" [" + synonym + "]");
			}
		}
	}

}
