package be.svlandeg.diffany.usecase.arabidopsis.osmotic;

import java.io.File;
import java.io.IOException;
import java.net.URI;
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
import be.svlandeg.diffany.usecase.ExpressionDataAnalysis;
import be.svlandeg.diffany.usecase.arabidopsis.CornetData;
import be.svlandeg.diffany.usecase.arabidopsis.GenePrinter;
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
	
	private URI ppi_file;
	private URI reg_file;
	
	/**
	 * The constructor defines a few properties of this analysis, such as where to fetch the PPI/regulatory data.
	 */
	public RunAnalysis()
	{
		// put to null when you don't want to consult these resources
		ppi_file = new CornetData().getCornetPPI();
		reg_file = new CornetData().getCornetReg();		
	}

	/**
	 * Run the full analysis pipeline.
	 * Currently, the data directory is hard coded to point to Sofie's D drive (TODO v2.1).
	 * 
	 * @param args these requirede CL arguments are currently not parsed
	 * @throws Exception whenever something goes wrong
	 */
	public static void main(String[] args) throws Exception
	{
		System.out.println("Performing osmotic data analysis");
		System.out.println("");

		//String inputRoot = "D:" + File.separator + "diffany-osmotic"; // Sofie @ PSB
		String inputRoot = "C:/Users/Sloffie/Documents/phd/diffany_data/osmotic"; // Sofie @ home

		File osmoticStressDir = new DataIO(inputRoot).getRootOsmoticStressDir();
		String outputDir = osmoticStressDir + File.separator + "output";
		
		boolean performStep1FromRaw = false;
		boolean performStep1FromSupplemental = true;
		boolean performStep2ToNetwork = true;
		boolean performStep3ToFile = false;
		boolean performStep4FromFile = false;
		
		if (performStep1FromRaw ==  performStep1FromSupplemental)
		{
			System.out.println("Select exactly one option to perform the first step of the analysis!");
			return;
		}

		double threshold_strict = 0.05;
		double threshold_fuzzy = 0.1;
		boolean writeHeaders = true;
		boolean allowVirtualEdges = false;
		
		boolean selfInteractions = false;
		boolean neighbours = true;
		//int min_neighbourcount = 1;		// TODO this is currently not implemented anymore
		boolean includeUnknownReg = false;
		
		String overexpressionFile = null;
		Set<InputNetwork> networks = null;
		RunAnalysis ra = new RunAnalysis();
		
		
		/* STEP 1 - OPTION 1: PROCESS RAW CELL DATA TO PRODUCE OVEREXPRESSION VALUES WITH DIFFANY */
		if (performStep1FromRaw)
		{
			System.out.println("1. Transforming CELL data into overexpression values");
			System.out.println("");
			String fileName = "differential_values.txt";
			overexpressionFile = ra.fromRawToOverexpression(osmoticStressDir, fileName);
			System.out.println(" written to " + overexpressionFile);
		}

		/* STEP 1 - OPTION 2: GET PUBLISHED OVEREXPRESSION VALUES FROM THE OSMOTIC PAPER */
		if (performStep1FromSupplemental)
		{
			System.out.println("1. Reading published overexpression values");
			System.out.println("");
			overexpressionFile = osmoticStressDir + File.separator + "clean_Inze_Supplemental_Dataset_1.tab";
		}

		/* STEP 2: USE OVEREXPRESSION VALUES TO CREATE NETWORKS */
		if (performStep2ToNetwork)
		{
			System.out.println("2. Transforming overexpression values into networks");
			System.out.println("   selfinteractions " + selfInteractions + " / neighbours " + neighbours + " / includeUnknownReg " + includeUnknownReg);
			System.out.println("");
			try
			{
				networks = ra.fromOverexpressionToNetworks(new File(overexpressionFile), 1, threshold_strict, threshold_fuzzy, selfInteractions, neighbours, includeUnknownReg);
			}
			catch (IllegalArgumentException e)
			{
				System.out.println("Trouble parsing from " + overexpressionFile + ": " + e.getMessage());
				return;
			}
			catch (IOException e)
			{
				System.out.println("IO trouble reading from " + overexpressionFile + ": " + e.getMessage());
				return;
			}
		}

		/* STEP 3: WRITE NETWORKS TO FILE */
		if (performStep2ToNetwork && performStep3ToFile)
		{
			System.out.println("");
			System.out.println("3. Writing output networks to " + outputDir);
			System.out.println("   allowVirtualEdges " + allowVirtualEdges);
			for (InputNetwork net : networks)
			{
				NetworkIO.writeNetworkToDir(net, net.getNodeMapper(), new File(outputDir, net.getName()), writeHeaders, allowVirtualEdges);
			}
		}

		/* STEP 4: READ NETWORKS BACK IN FROM FILE */
		if (performStep4FromFile)
		{
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
		}

		System.out.println("");
		System.out.println("Done!");
	}

	/**
	 * This first step in the pipeline processes raw .CELL files and produces a .tab file of the calculated p-values etc.
	 */
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
			deAnalysis.findDEGenes(exeR, outputFile);
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
	private Set<InputNetwork> fromOverexpressionToNetworks(File overExpressionFile, int firstID, double threshold_strict, double threshold_fuzzy, boolean selfInteractions, boolean neighbours, boolean includeUnknownReg) throws IOException, URISyntaxException
	{
		Set<InputNetwork> networks = new HashSet<InputNetwork>();
		GenePrinter gp = new GenePrinter();

		OverexpressionIO io = new OverexpressionIO();
		ExpressionDataAnalysis dataAn= new ExpressionDataAnalysis(gp);
		NetworkConstruction constr = new NetworkConstruction(gp);

		NodeMapper nm = new DefaultNodeMapper();
		EdgeOntology eo = new DefaultEdgeOntology();
		
		List<OverexpressionData> datasets = io.readDatasets(overExpressionFile, false);
		Set<String> all_nodeIDs_strict = new HashSet<String>();
		Set<String> all_nodeIDs_fuzzy = new HashSet<String>();
		
		for (OverexpressionData data : datasets)
		{
			System.out.println("");
			System.out.println(data.getName() + ": " + data.getArrayIDs().size() + " IDs analysed");

			Set<String> nodes_strict = nm.getNodeIDs(dataAn.getSignificantGenes(data, threshold_strict).keySet());
			all_nodeIDs_strict.addAll(nodes_strict);
			
			Set<String> nodes_fuzzy = nm.getNodeIDs(dataAn.getSignificantGenes(data, threshold_fuzzy).keySet());
			nodes_fuzzy.removeAll(nodes_strict);
			all_nodeIDs_fuzzy.addAll(nodes_fuzzy);
			
			System.out.println("  Found " + nodes_strict.size() + " differentially expressed genes at threshold " + threshold_strict + " and " + nodes_fuzzy.size() + " additional ones at threshold " + threshold_fuzzy);
		}
		
		// if they are strict DE once, they do not need to be in the fuzzy set also
		all_nodeIDs_fuzzy.removeAll(all_nodeIDs_strict);
		System.out.println("");
		System.out.println("Total: " + all_nodeIDs_strict.size() + " strict differentially expressed genes at threshold " + threshold_strict + " and " + all_nodeIDs_fuzzy.size() + " additional ones at threshold " + threshold_fuzzy);
	
		Set<String> expandedNetwork = constr.expandNetwork(nm, all_nodeIDs_strict, all_nodeIDs_fuzzy, ppi_file, reg_file, selfInteractions, neighbours, includeUnknownReg);
		System.out.println("  Found " + expandedNetwork.size() + " total nodes");
		
		Set<Node> all_nodes = new HashSet<Node>();
		for (String locusID : expandedNetwork)
		{
			String symbol = gp.getSymbolByLocusID(locusID);
			if (symbol == null)
			{
				symbol = locusID;
			}
			all_nodes.add(new Node(locusID, symbol, false));
		}
		
		
			
			/* TODO - refactor
			{
				Logger logger = new Logger();
				NetworkCleaning cleaning = new NetworkCleaning(logger);
				
				
				Set<Edge> edges = constr.createAllEdgesFromDiffData_old(nodes_strict, true, ppi_file, reg_file, selfInteractions, neighbours, includeUnknownReg);
				System.out.println("  Found " + edges.size() + " total edges");
	
				// The set of nodes is not updated at this point, but the Network constructor automatically adds the new ones from the edge information
				InputNetwork net = new InputNetwork(data.getName(), firstID++, new HashSet<Node>(nodes_strict.keySet()), edges, nm);
	
				System.out.println("  Cleaning network:");
	
				InputNetwork cleannet = cleaning.fullInputCleaning(net, nm, eo);
				networks.add(cleannet);
	
				for (LogEntry msg : logger.getAllLogMessages())
				{
					System.out.println(msg);
				}
	
				Set<Node> deNodes = nodes_strict.keySet();
				deNodes.retainAll(cleannet.getNodes());
				System.out.println("  Final network: " + cleannet.getEdges().size() + " non-redundant edges between " + cleannet.getNodes().size() + " nodes of which " + deNodes.size() + " DE nodes");
			}*/

		return networks;
	}
}
