package be.svlandeg.diffany.usecase.arabidopsis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

/**
 * This class is used to print all calculated overexpression values to one big file
 * 
 * @author Sofie Van Landeghem
 */
public class OverexpressionIO
{

	/**
	 * Method which will print all genes in the experiments and their values across all types of statistical comparisons
	 * @param filename the file location where all calculated p-values etc will be written
	 * @param datasets
	 * @param gp the gene printer which can fetch the locus IDs and synonyms for array IDs
	 * @throws IOException when the 
	 */
	public void printDatasets(String filename, List<OverexpressionData> datasets, GenePrinter gp) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));
		
		// HEADER
		writer.append("arrayID \t");
		for (OverexpressionData data : datasets)
		{
			String name = data.getName();
			writer.append("FDR" + name + "\t" + "FC" + name + "\t" + "p-value" + name + "\t");
		}
		writer.append("synonyms \t");
		writer.newLine();
		writer.flush();
		
		// DATA
		DecimalFormat df = new DecimalFormat("#.###");
		for (String arrayID: datasets.get(0).getArrayIDs())
		{
			writer.append(arrayID + " \t");
			for (OverexpressionData data : datasets)
			{
				double foldchange = data.getFoldchange(arrayID);
				double pvalue = data.getPvalue(arrayID);
				double FDR = data.getFDR(arrayID);
				
				writer.append(df.format(FDR) + " \t" + df.format(foldchange) + " \t" + df.format(pvalue) + " \t");
				writer.flush();
			}
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
