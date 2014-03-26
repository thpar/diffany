package be.svlandeg.diffany.usecase.arabidopsis.tf;

import java.io.File;
import java.io.IOException;


/**
 * This class provides the analysis pipeline and calls our procedures necessary
 * to perform the full analysis, from data input to result output.
 * 
 * @author Sofie Van Landeghem
 */
public class RunAnalysis
{
	
	/**
	 * Run the full analysis pipeline.
	 * Currently, the data directory is hard coded to point to Sofie's D drive (TODO v2.1).
	 * 
	 * @param args these requirede CL arguments are currently not parsed
	 */
	public static void main(String[] args)
	{
		System.out.println("Performing TF-target analysis");
		System.out.println("");
		String inputRoot = "D:" + File.separator + "diffany-tf";
		DataIO io = new DataIO(inputRoot);
		File tfTargetFile = io.getTFs();
		InputProcessing input = new InputProcessing();
		
		try
		{
			input.processTFData(tfTargetFile);
		}
		catch(IOException e)
		{
			String errorMsg = "Error reading input data from " + inputRoot + ": " + e.getMessage();
			System.out.println(errorMsg);
		}

		System.out.println("");
		System.out.println("Done!");
	}

}
