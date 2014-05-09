package be.svlandeg.diffany.usecase.arabidopsis.osmotic;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import be.svlandeg.diffany.r.ExecuteR;
import be.svlandeg.diffany.usecase.arabidopsis.GenePrinter;

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
	 * @param exeR the environment to execute R scripts
	 * @param osmoticStressDir the directory containing the experimental data
	 * 
	 * @throws URISyntaxException
	 * @throws IOException
	 */
	public void findDEGenes(ExecuteR exeR, File osmoticStressDir) throws URISyntaxException, IOException
	{
		// TODO V2.1: currently this assumes libs "affy", "limma", etc are pre-installed!
		URL script3URL = Thread.currentThread().getContextClassLoader().getResource("Rcode/FindDEgenes.R");
		System.out.println(" Executing script to find differentially expressed genes: " + script3URL);
		System.out.println(" (this may take a minute ... please be patient and do not interrupt the execution) ");
		exeR.executeScript(script3URL);
		System.out.println("");

		System.out.println(" Analysing data: ");

		String[] samples = exeR.getStringArray("samples");
		System.out.println("  Samples: " + samples.length);

		String[] probesets = exeR.getStringArray("probesets");
		System.out.println("  Probe sets: " + probesets.length);

		List<String> suffixes = new ArrayList<String>();
		suffixes.add("_stress");
		//suffixes.add("_stress_time_1.5");
		//suffixes.add("_stress_time_3");
		//suffixes.add("_stress_time_12");
		//suffixes.add("_stress_time_24");
		
		GenePrinter gp = new GenePrinter();
		
		printByResult(exeR, gp, suffixes, 50);
	}
	
	/**
	 * Private method which will print the top X genes per type of comparison
	 * 
	 * @param exeR the environment to execute R scripts
	 * @gp the gene printer which can fetch the locus IDs and synonyms for array IDs
	 * @param suffixes the suffixes used in the R script to subscript the "topIDs" objects
	 * @param printMax the number of results (array IDs) that should be printed (if there are that many)
	 */
	private void printByResult(ExecuteR exeR, GenePrinter gp, List<String> suffixes, int printMax)
	{
		System.out.println("");
		
		for (String suffix : suffixes)
		{			
			System.out.println("Results for: " + suffix);

			String[] topIDs = exeR.getStringArray("topIDs" + suffix);
			
			if (topIDs == null)
			{
				System.out.println("   error: NULL");
			}
			else
			{
				; //Integer.MAX_VALUE;
				int toPrint = Math.min(printMax, topIDs.length);

				System.out.println(" " + toPrint + " top most DE genes: ");
				System.out.println("");
				
				
				for (int i = 0; i < toPrint; i++)
				{
					String arrayID = topIDs[i];
					int rindex = i+1;
					
					double foldChange = exeR.getDoubleValue("toptable" + suffix + "[" + rindex + ",'logFC']");
					double pvalue = exeR.getDoubleValue("toptable" + suffix + "[" + rindex + ",'P.Value']");
					double FDR = exeR.getDoubleValue("toptable" + suffix + "[" + rindex + ",'adj.P.Val']");
					
					List<String> results = gp.getSynonyms(arrayID);
					System.out.println("  " + (i + 1) + "." + arrayID + " - FDR: " + FDR + " - FC: " + foldChange + " - pvalue: " + pvalue);
					for (String synonyms : results)
					{
						System.out.println("   \t" + synonyms);
					}
					System.out.println("");
				}
			}
			System.out.println(" ************************************************************** ");
		}
	}
}
