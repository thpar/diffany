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

		System.out.println(" Analysing data: ");

		String[] samples = exeR.getStringArray("samples");
		System.out.println("  Samples: " + samples.length);

		String[] probesets = exeR.getStringArray("probesets");
		System.out.println("  Probe sets: " + probesets.length);

		System.out.println("");
		GenePrinter gp = new GenePrinter();
		System.out.println("");

		List<String> suffixes = new ArrayList<String>();
		suffixes.add("_stress");
		//suffixes.add("_stress_time_1.5");
		//suffixes.add("_stress_time_3");
		//suffixes.add("_stress_time_12");
		//suffixes.add("_stress_time_24");

		for (String suffix : suffixes)
		{
			System.out.println("  results for: " + suffix);

			String[] topIDs = exeR.getStringArray("topIDs" + suffix);
			if (topIDs == null)
			{
				System.out.println("   error: NULL");
			}
			else
			{
				int printMax = 100; //Integer.MAX_VALUE;
				int toPrint = Math.min(printMax, topIDs.length);

				System.out.println(toPrint + " top most DE genes");

				for (int i = 0; i < toPrint; i++)
				{
					System.out.println("");
					String arrayID = topIDs[i];
					List<String> results = gp.printGene(arrayID);
					for (String result : results)
					{
						System.out.println("  " + (i + 1) + ".\t" + result);
					}
				}
			}
			System.out.println(" ******************************* ");
		}
	}
}
