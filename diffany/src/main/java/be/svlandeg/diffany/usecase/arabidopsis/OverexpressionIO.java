package be.svlandeg.diffany.usecase.arabidopsis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * This class is used to print all calculated overexpression values to one big file
 * 
 * @author Sofie Van Landeghem
 */
public class OverexpressionIO
{

	private static String arrayID_column = "arrayID";
	private static String FDR_column = "FDR";
	private static String FC_column = "FC";
	private static String pValue_column = "p-value";
	private static String synonyms_column = "synonyms";

	/**
	 * Print all genes in the experiments and their values across all types of statistical comparisons
	 * @param filename the file location where all calculated p-values etc will be written
	 * @param datasets the overexpression datasets
	 * @param gp the gene printer which can fetch the locus IDs and synonyms for array IDs
	 * @throws IOException when the file could not be written properly
	 */
	public void printDatasets(String filename, List<OverexpressionData> datasets, GenePrinter gp) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		// HEADER
		writer.append(arrayID_column + " \t");
		for (OverexpressionData data : datasets)
		{
			String name = data.getName();
			writer.append(FDR_column + name + "\t" + FC_column + name + "\t" + pValue_column + name + "\t");
		}
		writer.append(synonyms_column + " \t");
		writer.newLine();
		writer.flush();

		// DATA
		DecimalFormat df = new DecimalFormat("#.###");
		DecimalFormatSymbols symbols = df.getDecimalFormatSymbols();
		symbols.setDecimalSeparator(',');
		df.setDecimalFormatSymbols(symbols);
		
		boolean rawArrayIDs = datasets.get(0).indexedByRawArrayIDs();

		for (String ID : datasets.get(0).getArrayIDs())
		{
			writer.append(ID + " \t");
			for (OverexpressionData data : datasets)
			{
				double foldchange = data.getFoldchange(ID);
				double pvalue = data.getPvalue(ID);
				double FDR = data.getFDR(ID);

				writer.append(df.format(FDR) + " \t" + df.format(foldchange) + " \t" + df.format(pvalue) + " \t");
				writer.flush();
			}
			List<String> synonymList = gp.getSynonymsByArrayID(ID);
			if (!rawArrayIDs)
			{
				synonymList.clear();
				synonymList.add(gp.getSynonymsByLocusID(ID));
			}
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

	/**
	 * Read all genes and their corresponding p-values of the statistical comparisons defined by the header line
	 * 
	 * @param filename the file location where all calculated p-values etc are recorded
	 * @param rawArrayIDs whether or not this refers to raw array IDs (if not, locus tags are expected)
	 * @return the list of overexpression datasets
	 * @throws IOException when the file could not be read properly
	 * @throws IllegalArgumentException if the input file is not formatted as expected
	 */
	public List<OverexpressionData> readDatasets(String filename, boolean rawArrayIDs) throws IOException, IllegalArgumentException
	{
		return readDatasets(new File(filename), rawArrayIDs);
	}

	/**
	 * Read all genes and their corresponding p-values of the statistical comparisons defined by the header line
	 * @param inputfile the file which stores all calculated p-values etc 
	 * @param rawArrayIDs whether or not this refers to raw array IDs (if not, locus tags are expected)
	 * @return the list of overexpression datasets
	 * @throws IOException when the file could not be read properly
	 * @throws IllegalArgumentException if the input file is not formatted as expected
	 */
	public List<OverexpressionData> readDatasets(File inputfile, boolean rawArrayIDs) throws IOException, IllegalArgumentException
	{
		List<OverexpressionData> datasets = new ArrayList<OverexpressionData>();

		BufferedReader reader = new BufferedReader(new FileReader(inputfile));

		// READ AND CHECK HEADER
		String header = reader.readLine();
		//System.out.println("Read header: " + header);
		StringTokenizer header_stok = new StringTokenizer(header, "\t");

		boolean formatOK = true;

		String first_token = header_stok.nextToken().trim();
		if (!arrayID_column.equals(first_token))
		{
			formatOK = false;
		}
		String next_token = header_stok.nextToken().trim();
		while (formatOK && next_token != null && next_token.startsWith(FDR_column))
		{
			int splits = next_token.indexOf(FDR_column);
			String suffix = next_token.substring(splits + 3);

			next_token = header_stok.nextToken().trim();
			String expected_FC = FC_column + suffix;
			if (!expected_FC.equals(next_token))
			{
				formatOK = false;
			}

			next_token = header_stok.nextToken().trim();
			String expected_pValue = pValue_column + suffix;
			if (!expected_pValue.equals(next_token))
			{
				formatOK = false;
			}

			if (formatOK)
			{
				OverexpressionData data = new OverexpressionData(suffix, rawArrayIDs);
				datasets.add(data);
			}
			if (header_stok.hasMoreTokens())
			{
				next_token = header_stok.nextToken();
			}
			else
			{
				next_token = null;
			}
		}

		// TODO. Currently no check on synonyms column because it's not strictly necessary
		/*
		next_token = header_stok.nextToken();
		if (! synonyms_column.equals(next_token))
		{
			formatOK = false;
		}
		*/

		if (!formatOK || datasets.isEmpty())
		{
			reader.close();
			String errormsg = "Unexpected input file format!" + " Expected is a tab-delimited file with header line arrayID, followed by (FDR,FC,p-value) triples" + " for each comparison, and finally a synonyms column";
			throw new IllegalArgumentException(errormsg);
		}

		// READ ALL OTHER DATA
		String line = reader.readLine();
		//System.out.println("Read line: " + line);
		while (line != null)
		{
			// TODO v2.1: this is a quick hack to work around the european locale
			line = line.replace(',', '.');
			StringTokenizer stok = new StringTokenizer(line, "\t");
			String arrayID = stok.nextToken().trim();

			// Read all triples for this array ID, i.e. 3 consecutive columns define 1 comparison
			for (OverexpressionData data : datasets)
			{
				double FDR = Double.parseDouble(stok.nextToken().trim());
				double FC = Double.parseDouble(stok.nextToken().trim());
				double pvalue = Double.parseDouble(stok.nextToken().trim());
				data.addResult(arrayID, FC, pvalue, FDR);
			}
			line = reader.readLine();
		}

		reader.close();
		return datasets;
	}

}
