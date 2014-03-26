package be.svlandeg.diffany.usecase.arabidopsis.tf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

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
		System.out.println(" Reading " + path + ":");

		BufferedReader reader = new BufferedReader(new FileReader(TFtargetFile));
		String line = reader.readLine();
		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			String TF = stok.nextToken();
			String target = stok.nextToken();

			System.out.println(TF + " --> " + target);

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
	public void processExpressionData(File expDir) throws IOException
	{
		for (File expFile : expDir.listFiles())
		{
			String path = expFile.getAbsolutePath();
			System.out.println(" Reading " + path + ":");

			BufferedReader reader = new BufferedReader(new FileReader(expFile));
			String line = reader.readLine();
			while (line != null)
			{
				StringTokenizer stok = new StringTokenizer(line, "\t");
				String gene = stok.nextToken();
				List<String> values = new ArrayList<String>();
				while (stok.hasMoreTokens())
				{
					values.add(stok.nextToken());
				}
				System.out.println(gene + " --> " + values.size() + " values");

				line = reader.readLine();
			}
			reader.close();
			System.out.println(" ");
		}
	}

}
