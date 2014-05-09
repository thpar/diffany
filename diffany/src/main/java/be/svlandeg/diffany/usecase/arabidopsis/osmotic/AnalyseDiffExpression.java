package be.svlandeg.diffany.usecase.arabidopsis.osmotic;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

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
	 * @throws IOException when an IO error occurs
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
		suffixes.add("_mannitol");
		suffixes.add("_1.5");
		suffixes.add("_3");
		suffixes.add("_12");
		suffixes.add("_24");

		GenePrinter gp = new GenePrinter();

		printByStatTest(exeR, gp, suffixes, 10);

		String outputValues = osmoticStressDir + File.separator + "differential_values.txt";
		printByGene(exeR, gp, suffixes, outputValues);
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

					List<String> synonymList = gp.getSynonyms(arrayID);
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
	 * Private method which will print all genes in the experiments and their values across all types of statistical comparisons.
	 * 
	 * @param exeR the environment to execute R scripts
	 * @@param gp the gene printer which can fetch the locus IDs and synonyms for array IDs
	 * @param suffixes the suffixes used in the R script to subscript the "topIDs" objects
	 * @param outputValues the file location where all calculated p-values etc will be written
	 * @throws IOException when an IO error occurs
	 */
	private void printByGene(ExecuteR exeR, GenePrinter gp, List<String> suffixes, String outputValues) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputValues));
		System.out.println("");

		SortedSet<String> arrayIDs = new TreeSet<String>();
		Map<String, Double> foldchanges = new HashMap<String, Double>();
		Map<String, Double> pvalues = new HashMap<String, Double>();
		Map<String, Double> FDRs = new HashMap<String, Double>();

		System.out.println("All results are written to " + outputValues);

		writer.append("arrayID \t");

		for (String suffix : suffixes)
		{
			writer.append("FDR" + suffix + "\t" + "FC" + suffix + "\t" + "p-value" + suffix + "\t");
			String[] topIDs = exeR.getStringArray("topIDs" + suffix);

			if (topIDs != null)
			{
				for (int i = 0; i < topIDs.length; i++)
				{
					String arrayID = topIDs[i];
					arrayIDs.add(arrayID);

					String key = suffix + "*" + arrayID;

					if (foldchanges.containsKey(key) || foldchanges.containsKey(key) || foldchanges.containsKey(key))
					{
						System.out.println(" WARNING: encountered " + arrayIDs + " twice in " + suffix + " ?! ");
					}

					int rindex = i + 1;

					foldchanges.put(key, exeR.getDoubleValue("toptable" + suffix + "[" + rindex + ",'logFC']"));
					pvalues.put(key, exeR.getDoubleValue("toptable" + suffix + "[" + rindex + ",'P.Value']"));
					FDRs.put(key, exeR.getDoubleValue("toptable" + suffix + "[" + rindex + ",'adj.P.Val']"));
				}
			}
		}
		writer.append("synonyms \t");
		writer.newLine();
		writer.flush();

		DecimalFormat df = new DecimalFormat("#.###");

		for (String arrayID : arrayIDs)
		{
			writer.append(arrayID+ "\t");
			for (String suffix : suffixes)
			{
				String key = suffix + "*" + arrayID;
				writer.append(df.format(FDRs.get(key)) + "\t" + df.format(foldchanges.get(key)) + "\t" + df.format(foldchanges.get(key)) + "\t");
			}
			writer.flush();

			List<String> synonymList = gp.getSynonyms(arrayID);

			for (String synonyms : synonymList)
			{
				writer.append(synonyms + " /// ");
			}
			writer.newLine();
			writer.flush();
		}
		writer.flush();
		writer.close();
	}
}
