package be.svlandeg.diffany.usecase.arabidopsis.tf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;

import be.svlandeg.diffany.core.expression.ExpressionData;

/**
 * This class reads and processes the raw input TF-target data.
 * 
 * @author Sofie Van Landeghem
 */
public class InputProcessing
{

	/**
	 * Process the TF-target data
	 * 
	 * @param TFtargetFile the file containing the TF-target interactions
	 * @throws IOException
	 */
	public void processTFData(File TFtargetFile) throws IOException
	{
		String path = TFtargetFile.getAbsolutePath();
		System.out.println(" Reading " + path );

		BufferedReader reader = new BufferedReader(new FileReader(TFtargetFile));
		String line = reader.readLine();
		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			String TF = stok.nextToken();
			String target = stok.nextToken();

			//System.out.println(TF + " --> " + target);

			line = reader.readLine();
		}
		reader.close();
	}

	/**
	 * Process the expression data
	 * 
	 * @param expDir the directory containing the pre-processed expression files
	 * @throws IOException
	 */
	public Set<ExpressionData> processExpressionData(File expDir) throws IOException
	{
	    Set<ExpressionData> datasets = new HashSet<ExpressionData>();
	    
		for (File expFile : expDir.listFiles())
		{
			String path = expFile.getAbsolutePath();
			System.out.println(" Reading " + path);
			
			List<String> samples = new ArrayList<String>();
			List<String> genes = new ArrayList<String>();
			List<List<Double>> matrix = new ArrayList<List<Double>>();
			
			BufferedReader reader = new BufferedReader(new FileReader(expFile));
			
			String header = reader.readLine();
			StringTokenizer header_stok = new StringTokenizer(header, "\t");
			
			String first = header_stok.nextToken().trim();
			if (!first.equals("#"))
			{
				reader.close();
				String errormsg = "Unexpected format for the MA data, header line: " + header;
				throw new IllegalArgumentException(errormsg);
			}
			
			while (header_stok.hasMoreTokens())
			{
				samples.add(header_stok.nextToken());
			}
			
			String line = reader.readLine();
			while (line != null)
			{
				StringTokenizer stok = new StringTokenizer(line, "\t");
				String gene = stok.nextToken();
				genes.add(gene);
				
				List<Double> values = new ArrayList<Double>();
				while (stok.hasMoreTokens())
				{
					values.add(Double.parseDouble(stok.nextToken()));
				}
				
				matrix.add(values);
				line = reader.readLine();
				
			}
			
			int nr_genes = genes.size();
			int nr_samples = samples.size();
				
			double[][] expvalues = new double[nr_genes][nr_samples];
			for (int i = 0; i < nr_genes; i++)
			{
				List<Double> values = matrix.get(i);
				for (int j = 0; j < nr_samples; j++)
				{
					expvalues[i][j] = values.get(j);
				}
			}
 			
			ExpressionData dataset = new ExpressionData(expFile.getName(), genes, samples, expvalues, true);
			datasets.add(dataset);
			
			reader.close();
			System.out.println(" ");
		}
		return datasets;
	}
	

}
