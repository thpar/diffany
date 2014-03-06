package be.svlandeg.diffany.usecase.osmotic;

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
	 * Currently, the data directory is hard coded to point to Sofie's D drive (TODO).
	 * 
	 * @param args these requirede CL arguments are currently not parsed
	 */
	public static void main(String[] args)
	{
		String inputRoot = "D:" + File.separator + "diffany-osmotic";
		DataIO io = new DataIO(inputRoot);
		InputProcessing input = new InputProcessing();
		
		File osmoticStressDir = io.getOsmoticStressDir();
		try
		{
			input.processOsmoticData(osmoticStressDir);
		}
		catch(IOException e)
		{
			String errorMsg = "Error reading input data from " + inputRoot + ": " + e.getMessage();
			System.out.println(errorMsg);
		}
	}

}
