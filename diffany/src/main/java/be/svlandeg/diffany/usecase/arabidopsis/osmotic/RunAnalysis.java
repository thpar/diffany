package be.svlandeg.diffany.usecase.arabidopsis.osmotic;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.algorithms.NetworkCleaning;
import be.svlandeg.diffany.core.io.NetworkIO;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.project.LogEntry;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.core.semantics.EdgeOntology;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.r.ExecuteR;
import be.svlandeg.diffany.r.RBridge;
import be.svlandeg.diffany.usecase.arabidopsis.NetworkConstruction;
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
	 * @throws URISyntaxException
	 */
	public static void main(String[] args) throws IOException, URISyntaxException
	{
		System.out.println("Performing osmotic data analysis");
		System.out.println("");

		//String inputRoot = "D:" + File.separator + "diffany-osmotic";					// Sofie @ PSB
		String inputRoot = "C:/Users/Sloffie/Documents/phd/diffany_data/osmotic"; // Sofie @ home

		File osmoticStressDir = new DataIO(inputRoot).getRootOsmoticStressDir();
		RunAnalysis ra = new RunAnalysis();

		/*
		 * STEP 1 - OPTION 1: PROCESS RAW CELL DATA TO PRODUCE OVEREXPRESSION VALUES WITH DIFFANY
		 * 
		System.out.println("1. Transforming CELL data into overexpression values");
		System.out.println("");
		String fileName = "differential_values.txt";
		String overexpressionFile = ra.fromRawToOverexpression(osmoticStressDir, fileName);
		*/

		/*
		 * STEP 1 - OPTION 2: GET PUBLISHED OVEREXPRESSION VALUES FROM THE OSMOTIC PAPER
		 */
		System.out.println("1. Reading published overexpression values");
		System.out.println("");

		String overexpressionFile = osmoticStressDir + File.separator + "clean_Inze_Supplemental_Dataset_1.tab";

		/*
		 * STEP 2: USE OVEREXPRESSION VALUES TO CREATE NETWORKS
		 */
		System.out.println("2. Transforming overexpression values into networks");
		
		boolean selfInteractions = false;
		boolean neighbours = true;
		int min_neighbourcount = 1;
		boolean addPPI = true;
		boolean addReg = true;
		
		System.out.println("   selfinteractions " + selfInteractions + " / neighbours " + neighbours + " / min_neighbourcount " + min_neighbourcount + " / addPPI " + addPPI + " / addReg " + addReg);
		System.out.println("");
		
		Set<InputNetwork> networks = ra.fromOverexpressionToNetworks(new File(overexpressionFile), 0.05, selfInteractions, neighbours, min_neighbourcount, addPPI, addReg);

		/*
		 * STEP 3: WRITE NETWORKS TO FILE
		 */
		String outputDir = osmoticStressDir + File.separator + "output";
		boolean writeHeaders = true; 
		boolean allowVirtualEdges = false;

		System.out.println("");
		System.out.println("3. Writing output networks to " + outputDir);
		System.out.println("   allowVirtualEdges " + allowVirtualEdges);

		for (InputNetwork net : networks)
		{
			NetworkIO.writeNetworkToDir(net, net.getNodeMapper(), new File(outputDir, net.getName()), writeHeaders, allowVirtualEdges);
		}

		/*
		 * STEP 4: READ NETWORKS BACK IN FROM FILE
		 */
		System.out.println("");
		System.out.println("4. Reading networks from " + outputDir);

		Set<InputNetwork> readNetworks = NetworkIO.readGenericInputNetworksFromSubdirs(new File(outputDir), new DefaultNodeMapper(), writeHeaders);

		
		for (InputNetwork rn : readNetworks)
		{
			System.out.println("");
			System.out.println(" " + rn.getStringRepresentation() + ": ");
			System.out.print(" " + rn.getNodes().size() + " nodes and " + rn.getEdges().size() + " edges");
			System.out.println(" (" + rn.getNodesByVirtualState(true).size() + " virtual nodes)");
		}

		System.out.println("");
		System.out.println("Done!");
	}

	/**
	 * This first step in the pipeline processes raw .CELL files and produces a .tab file of the calculated p-values etc.
	 */
	@SuppressWarnings("unused")
	private String fromRawToOverexpression(File osmoticStressDir, String overExpressionFile)
	{
		InputProcessing input = new InputProcessing();
		AnalyseDiffExpression deAnalysis = new AnalyseDiffExpression();

		String outputLog = osmoticStressDir + File.separator + "R_log.txt"; // can also be null
		String outputFile = osmoticStressDir + File.separator + overExpressionFile;

		RBridge bridge = new RBridge(outputLog);
		try
		{
			ExecuteR exeR = new ExecuteR(bridge);
			input.processOsmoticCELLData(exeR, osmoticStressDir);
			deAnalysis.findDEGenes(exeR, osmoticStressDir, outputFile);
		}
		catch (IOException e)
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

		if (!errors.isEmpty())
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

	/**
	 * Second step in the pipeline: use the overexpression values to generate networks
	 * 
	 * @throws URISyntaxException
	 */
	private Set<InputNetwork> fromOverexpressionToNetworks(File overExpressionFile, double threshold, boolean selfInteractions, boolean neighbours, int min_neighbourcount, boolean addPPI, boolean addReg) throws IOException, URISyntaxException
	{
		Set<InputNetwork> networks = new HashSet<InputNetwork>();
		
		OverexpressionIO io = new OverexpressionIO();
		NetworkConstruction constr = new NetworkConstruction();
		
		NodeMapper nm = new DefaultNodeMapper();
		EdgeOntology eo = new DefaultEdgeOntology();

		List<OverexpressionData> datasets = io.readDatasets(overExpressionFile, false);
		for (OverexpressionData data : datasets)
		{
			Logger logger = new Logger();
			NetworkCleaning cleaning = new NetworkCleaning(logger);
			
			System.out.println("");
			System.out.println(data.getName() + ": " + data.getArrayIDs().size() + " IDs analysed");

			Map<Node, Double> nodes = constr.getSignificantGenes(data, threshold);
			System.out.println("  Found " + nodes.size() + " differentially expressed genes");

			Set<Edge> edges = constr.constructVirtualRegulations(nodes);
			System.out.println("  Found " + edges.size() + " virtual regulations");

			if (addPPI)
			{
				Set<Edge> PPIedges = constr.readPPIsByLocustags(nodes.keySet(), selfInteractions, neighbours, min_neighbourcount);
				System.out.println("  Found " + PPIedges.size() + " PPIs");
				edges.addAll(PPIedges);
			}
			
			if (addPPI)
			{
				Set<Edge> regEdges = constr.readPPIsByLocustags(nodes.keySet(), selfInteractions, neighbours, min_neighbourcount);
				System.out.println("  Found " + regEdges.size() + " regulations");
				edges.addAll(regEdges);
			}

			System.out.println("  Found " + edges.size() + " total edges");

			InputNetwork net = new InputNetwork(data.getName(), new HashSet<Node>(nodes.keySet()), edges, nm);
			
			System.out.println("  Cleaning network:");

			InputNetwork cleannet = cleaning.fullInputCleaning(net, nm, eo);
			networks.add(cleannet);
			
			for (LogEntry msg : logger.getAllLogMessages())
			{
				System.out.println(msg);
			}
			
			System.out.println("  Final network: " + cleannet.getEdges().size() + " non-redundant edges of which " + cleannet.getEdgesByVirtualState(true).size() + " virtual edges");
		}
		
		
		
		return networks;
	}
}
