package be.svlandeg.diffany.usecase.arabidopsis.osmotic;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.r.ExecuteR;
import be.svlandeg.diffany.usecase.arabidopsis.MapID;

/**
 * This class processes normalized expression data, identifying differentially expressed genes
 * and differential coexpression links, with R scripts.
 * 
 * @author Sofie Van Landeghem
 */
public class AnalyseDiffExpression
{
	
	/**
	 * Identify differentially expressed (DE) genes, with R scripts.
	 * 
	 * Currently, the needed R script is loaded from the context, and is defined in the 'resources' folder of the Maven project.
	 * TODO v2.1: will this code work when packaged inside a jar or will we need to create a tmp file?
	 * 
	 * @throws URISyntaxException
	 * @throws IOException
	 */
	public void findDEGenes(ExecuteR exeR, File osmoticStressDir) throws URISyntaxException, IOException
	{
		// TODO V2.1: currently this assumes libs "affy", "affyPLM" and "org.Dm.eg.db" are pre-installed!
		URL script3URL = Thread.currentThread().getContextClassLoader().getResource("Rcode/FindDEgenes.R");
		System.out.println(" Executing script to find differentially expressed genes: " + script3URL);
		System.out.println(" (this may take a minute ... please be patient and do not interrupt the execution) ");
		exeR.executeScript(script3URL);
		System.out.println("");
		
		URL arrayMappingURL = Thread.currentThread().getContextClassLoader().getResource("data/affy_ATH1_ID_mapping.tab");
		System.out.println(" Fetching array ID mapping data: " + arrayMappingURL);
		Map<String, Set<String>> arrayidmapping = new MapID().getAllArrayMappings(new File(arrayMappingURL.toURI())); 
		
		URL locusMappingURL = Thread.currentThread().getContextClassLoader().getResource("data/TAIR10_NCBI_GENEID_mapping.tab");
		System.out.println(" Fetching Locus ID mapping data: " + locusMappingURL);
		Map<String, String> locusidmapping = new MapID().getLocusGIDMappings(new File(locusMappingURL.toURI()));
		
		URL symbolMappingURL = Thread.currentThread().getContextClassLoader().getResource("data/EVEX_synonyms_3702.tab");
		System.out.println(" Fetching gene symbol mapping data: " + symbolMappingURL);
		Map<String, String> symbolmapping = new MapID().getSymbolMappings(new File(symbolMappingURL.toURI()));
		
		System.out.println("");
		System.out.println(" Analysing data: ");
		
		String[] samples = exeR.getStringArray("samples");
		System.out.println("  Samples: " + samples.length);
		
		String[] probesets = exeR.getStringArray("probesets");
		System.out.println("  Probe sets: " + probesets.length);
		
		//double[][] expressionValues = exeR.getDoubleMatrix("expressionMatrix");
		//System.out.println("  Expression values dimension: " + expressionValues.length + " - " + expressionValues[0].length);
		
		System.out.println("");
		String[] topIDs = exeR.getStringArray("topIDs");
		System.out.println(" Top most DE genes: ");
		for (int i = 0; i < topIDs.length; i++)
		{
			String arrayID = topIDs[i];
			Set<String> locusIDs = arrayidmapping.get(arrayID);
			for (String locusID : locusIDs)
			{
				String egid = locusidmapping.get(locusID);
				String symbol = symbolmapping.get(egid);
				System.out.print("  " + (i+1) + ". " + arrayID + " - " + locusID + " - GID:" + egid);
				if (symbol != null && ! symbol.equals(locusID))
				{
					System.out.print(" - " + symbol);
				}
				System.out.println("");
			}
		}
		System.out.println("");
	}

}
