package be.svlandeg.diffany.usecase.arabidopsis.osmotic;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import be.svlandeg.diffany.r.ExecuteR;
import be.svlandeg.diffany.usecase.arabidopsis.GenePrinter;
import be.svlandeg.diffany.usecase.arabidopsis.OverexpressionData;
import be.svlandeg.diffany.usecase.arabidopsis.OverexpressionIO;

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
	 * This method can write an output file containing all p-values etc, which will be called 'differential_values.txt' and placed within the 
	 * 'osmoticStressDir' root directory.
	 * 
	 * Currently, the needed R script is loaded from the context, and is defined in the 'resources' folder of the Maven project.
	 * TODO v2.1: will this code work when packaged inside a jar or will we need to create a tmp file?
	 * 
	 * @param exeR the environment to execute R scripts
	 * @param osmoticStressDir the directory containing the experimental data
	 * @param outputValues the location to which to write the results (if null, no results will be written)
	 * 
	 * @throws URISyntaxException when the R code can not be read properly
	 * @throws IOException when an IO error occurs
	 */
	public void findDEGenes(ExecuteR exeR, File osmoticStressDir, String outputValues) throws URISyntaxException, IOException
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
		//suffixes.add("_mannitol");
		suffixes.add("_1.5");
		suffixes.add("_3");
		suffixes.add("_12");
		suffixes.add("_24");

		GenePrinter gp = new GenePrinter();

		printByStatTest(exeR, gp, suffixes, 3);

		if (outputValues != null)
		{
			System.out.println("");
			System.out.println("All results are being written to " + outputValues);
			System.out.println(" This may take a few minutes, but please do not interfere with the file during this process");
			
			List<OverexpressionData> datasets = assembleDataset(exeR, suffixes);
			OverexpressionIO io = new OverexpressionIO();
			io.printDatasets(outputValues, datasets, gp);
		}
	}

	/**
	 * Private method which will print the top X genes per type of comparison
	 * 
	 * @param exeR the environment to execute R scripts
	 * @gp the gene printer which can fetch the locus IDs and synonyms for array IDs
	 * @param suffixes the suffixes used in the R script to subscript the "topIDs" objects
	 * @param printMax the number of results (array IDs) that should be printed (if there are that many)
	 */
	private void printByStatTest(ExecuteR exeR, GenePrinter gp, List<String> suffixes, int printMax)
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
				int toPrint = Math.min(printMax, topIDs.length);

				System.out.println(" " + toPrint + " top most DE genes: ");
				System.out.println("");

				for (int i = 0; i < toPrint; i++)
				{
					String arrayID = topIDs[i];
					int rindex = i + 1;

					double foldChange = exeR.getDoubleValue("toptable" + suffix + "[" + rindex + ",'logFC']");
					double pvalue = exeR.getDoubleValue("toptable" + suffix + "[" + rindex + ",'P.Value']");
					double FDR = exeR.getDoubleValue("toptable" + suffix + "[" + rindex + ",'adj.P.Val']");

					List<String> synonymList = gp.getSynonymsByArrayID(arrayID);
					System.out.println("  " + (i + 1) + "." + arrayID + " - FDR: " + FDR + " - FC: " + foldChange + " - pvalue: " + pvalue);
					for (String synonyms : synonymList)
					{
						System.out.println("   \t" + synonyms);
					}
					System.out.println("");
				}
			}
			System.out.println(" ************************************************************** ");
		}
	}

	/**
	 * Private method which will assemble all overexpression values across all types of statistical comparisons, ready for printing.
	 * 
	 * @param exeR the environment to execute R scripts
	 * @param suffixes the suffixes used in the R script to subscript the "topIDs" objects
	 */
	private List<OverexpressionData> assembleDataset(ExecuteR exeR, List<String> suffixes)
	{
		List<OverexpressionData> datasets = new ArrayList<OverexpressionData>();
		for (String suffix : suffixes)
		{
			OverexpressionData data = new OverexpressionData(suffix, true);
			
			String[] topIDs = exeR.getStringArray("topIDs" + suffix);

			if (topIDs != null)
			{
				for (int i = 0; i < topIDs.length; i++)
				{
					String arrayID = topIDs[i];

					int rindex = i + 1;

					double foldchange = exeR.getDoubleValue("toptable" + suffix + "[" + rindex + ",'logFC']");
					double pvalue = exeR.getDoubleValue("toptable" + suffix + "[" + rindex + ",'P.Value']");
					double FDR = exeR.getDoubleValue("toptable" + suffix + "[" + rindex + ",'adj.P.Val']");

					data.addResult(arrayID, foldchange, pvalue, FDR);
				}
			}
			datasets.add(data);
		}
		return datasets;
		
	}
}
