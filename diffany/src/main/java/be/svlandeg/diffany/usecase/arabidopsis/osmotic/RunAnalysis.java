package be.svlandeg.diffany.usecase.arabidopsis.osmotic;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.List;

import be.svlandeg.diffany.r.ExecuteR;
import be.svlandeg.diffany.r.RBridge;
import be.svlandeg.diffany.usecase.arabidopsis.OverexpressionData;
import be.svlandeg.diffany.usecase.arabidopsis.OverexpressionIO;


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
	public static void main(String[] args) throws IOException
	{
		System.out.println("Performing osmotic data analysis");
		System.out.println("");
		
		String inputRoot = "D:" + File.separator + "diffany-osmotic";					// Sofie @ PSB
		//String inputRoot = "C:/Users/Sloffie/Documents/phd/diffany_data/osmotic";		// Sofie @ home
		
		DataIO dataIO = new DataIO(inputRoot);
		OverexpressionIO io = new OverexpressionIO();
		
		File osmoticStressDir = dataIO.getRootOsmoticStressDir();
		
		RunAnalysis ra = new RunAnalysis();
		
		/*
		 * OPTION 1: PROCESS RAW CELL DATA TO PRODUCE OVEREXPRESSION VALUES WITH DIFFANY
		 * 
		System.out.println("1. Transforming CELL data into overexpression values");
		System.out.println("");
		String fileName = "differential_values.txt";
		String overexpressionFile = ra.fromRawToOverexpression(osmoticStressDir, fileName);
		*/
		
		/*
		 * OPTION 2: USE PUBLISHED OVEREXPRESSION VALUES FROM THE OSMOTIC PAPER
		 */
		System.out.println("1. Reading published overexpression values");
		System.out.println("");
		
		String overexpressionFile = osmoticStressDir + File.separator + "clean_Inze_Supplemental_Dataset_1.tab";
		
		System.out.println("2. Transforming overexpression values into networks");
		System.out.println("");
		
		List<OverexpressionData> datasets = io.readDatasets(overexpressionFile);
		for (OverexpressionData data : datasets)
		{
			System.out.println(data.getName() + ": " + data.getArrayIDs().size() + " IDs");
		}
		
		System.out.println("");
		System.out.println("Done!");
	}
	
	/**
	 * This (private) step in the pipeline processes raw .CELL files and produces a .tab file of the calculated p-values etc.
	 */
	private String fromRawToOverexpression(File osmoticStressDir, String fileName)
	{
		InputProcessing input = new InputProcessing();
		AnalyseDiffExpression deAnalysis = new AnalyseDiffExpression();
		
		
		String outputLog = osmoticStressDir + File.separator + "R_log.txt";	// can also be null
		String outputFile = osmoticStressDir + File.separator + fileName;
		
		RBridge bridge = new RBridge(outputLog);
		try
		{ 
			ExecuteR exeR = new ExecuteR(bridge);
			input.processOsmoticCELLData(exeR, osmoticStressDir);
			deAnalysis.findDEGenes(exeR, osmoticStressDir, outputFile);
		}
		catch(IOException e)
		{
			String errorMsg = "Error reading input data from " + osmoticStressDir + ": " + e.getMessage();
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
		return outputFile;
	}
}
