package be.svlandeg.diffany.usecase.arabidopsis.osmotic;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.List;

import be.svlandeg.diffany.r.ExecuteR;
import be.svlandeg.diffany.r.RBridge;


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
		System.out.println("Performing osmotic data analysis");
		System.out.println("");
		String inputRoot = "D:" + File.separator + "diffany-osmotic";					// Sofie @ PSB
		//String inputRoot = "C:/Users/Sloffie/Documents/phd/diffany_data/osmotic";		// Sofie @ home 
		DataIO io = new DataIO(inputRoot);
		InputProcessing input = new InputProcessing();
		AnalyseDiffExpression deAnalysis = new AnalyseDiffExpression();
		
		File osmoticStressDir = io.getRootOsmoticStressDir();
		String outputLog = osmoticStressDir + File.separator + "R_log.txt";	// can also be null
		boolean printDetails = true;
		
		RBridge bridge = new RBridge(outputLog);
		try
		{ 
			ExecuteR exeR = new ExecuteR(bridge);
			input.processOsmoticCELLData(exeR, osmoticStressDir);
			deAnalysis.findDEGenes(exeR, osmoticStressDir, printDetails);
		}
		catch(IOException e)
		{
			String errorMsg = "Error reading input data from " + inputRoot + ": " + e.getMessage();
			System.out.println(errorMsg);
		}
		catch (URISyntaxException e)
		{
			System.out.println("Couldn't read R script : " + e.getMessage());
		}
		System.out.println("");
		List<String> errors = bridge.getErrorsFromLogfile();
		
		if (! errors.isEmpty())
		{
			System.out.println(" Errors occurred during R execution:");
			for (String error : errors)
			{
				System.out.println("  ! " + error);
			}
		}
		bridge.close();
		
		System.out.println("");
		System.out.println("Done!");
	}
}
