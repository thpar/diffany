package be.svlandeg.diffany.usecase.arabidopsis.osmotic;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.algorithms.NetworkCleaning;
import be.svlandeg.diffany.core.io.NetworkIO;
import be.svlandeg.diffany.core.networks.Condition;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.core.semantics.EdgeOntology;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;
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
		System.out.println("Performing osmotic data analysis - " + new Date());
		System.out.println("");

		//String inputRoot = "D:" + File.separator + "diffany-osmotic"; // Sofie @ PSB
		String inputRoot = "C:/Users/Sloffie/Documents/phd/diffany_data/osmotic"; // Sofie @ home

		File osmoticStressDir = new DataIO(inputRoot).getRootOsmoticStressDir();
		String outputDir = osmoticStressDir + File.separator + "output";
		
		boolean performStep1FromRaw = false;
		boolean performStep1FromSupplemental = false;
		boolean performStep2ToNetwork = false;
		boolean performStep3ToFile = false;
		boolean performStep4FromFile = true;
		boolean performStep5OneagainstAll = true;
		
		if (performStep1FromRaw == performStep1FromSupplemental && performStep2ToNetwork)
		{
			System.out.println("Select exactly one option to perform the first step of the analysis!");
			return;
		}
		if (performStep3ToFile && ! performStep2ToNetwork)
		{
			System.out.println("Can not perform the third step without the second!");
			return;
		}
		if (performStep5OneagainstAll && ! performStep4FromFile)
		{
			System.out.println("Can not perform the fourth step without the third!");
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
		
		double weight_cutoff = 0;
		
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

		Set<ConditionNetwork> conditionNets = new HashSet<ConditionNetwork>();
		ReferenceNetwork refNet = null;
		
		/* STEP 4: READ NETWORKS BACK IN FROM FILE */
		if (performStep4FromFile)
		{
			System.out.println("");
			System.out.println("4. Reading networks from " + outputDir);
			System.out.println("");
	
			Set<InputNetwork> readNetworks = NetworkIO.readGenericInputNetworksFromSubdirs(new File(outputDir), new DefaultNodeMapper(), writeHeaders);
			for (InputNetwork net : readNetworks)
			{
				if (net instanceof ReferenceNetwork)
				{
					if (refNet != null)
					{
						System.out.println(" Found more than 1 reference network ?! ");
						return;
					}
					refNet = (ReferenceNetwork) net;
				}
				else if (net instanceof ConditionNetwork)
				{
					conditionNets.add((ConditionNetwork) net);
				}
				else
				{
					System.out.println(" Found a strange input network: " + net);
				}
			}
			System.out.println(" Found a reference network: ");
			System.out.print(" " + refNet.getStringRepresentation() + ": ");
			System.out.println(" " + refNet.getNodes().size() + " nodes and " + refNet.getEdges().size() + " edges");
			System.out.println("");
			
			System.out.println(" Found " + conditionNets.size() + " condition-specific networks:");
			for (ConditionNetwork cn : conditionNets)
			{
				System.out.print(" " + cn.getStringRepresentation() + ": ");
				System.out.println(" " + cn.getNodes().size() + " nodes and " + cn.getEdges().size() + " edges");
			}
		}
		if (performStep5OneagainstAll)
		{
			System.out.println("");
			System.out.println("5. Performing one against all analysis at cutoff " + weight_cutoff);
			System.out.println("");
			if (refNet == null || conditionNets == null || conditionNets.isEmpty())
			{
				System.out.println(" Did not find the correct reference and condition-specific networks! ");
				return;
			}
			ra.runDiffany(refNet, conditionNets, weight_cutoff);
		}

		System.out.println("");
		System.out.println("Done! - " + new Date());
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
	 * Second step in the pipeline: use the overexpression values to generate networks: 
	 * 1 reference network + 1 condition-dependent network for each overexpression dataset that is read from the input
	 */
	private Set<InputNetwork> fromOverexpressionToNetworks(File overExpressionFile, int firstID, double threshold_strict, double threshold_fuzzy, boolean selfInteractions, boolean neighbours, boolean includeUnknownReg) throws IOException, URISyntaxException
	{
		Set<InputNetwork> networks = new HashSet<InputNetwork>();
		GenePrinter gp = new GenePrinter();

		OverexpressionIO io = new OverexpressionIO();
		ExpressionDataAnalysis dataAn= new ExpressionDataAnalysis();
		NetworkConstruction constr = new NetworkConstruction(gp);

		NodeMapper nm = new DefaultNodeMapper();
		EdgeOntology eo = new DefaultEdgeOntology();
		Logger logger = new Logger();
		NetworkCleaning cleaning = new NetworkCleaning(logger);
		
		List<OverexpressionData> datasets = io.readDatasets(overExpressionFile, false);
		
		Set<String> all_nodeIDs_strict = new HashSet<String>();
		Set<String> all_nodeIDs_fuzzy = new HashSet<String>();
		
		// Read all different experiments and determine their overexpressed genes
		for (OverexpressionData data : datasets)
		{
			System.out.println("");
			System.out.println(data.getName() + ": " + data.getArrayIDs().size() + " IDs analysed");

			Set<String> nodes_strict = dataAn.getSignificantGenes(data, threshold_strict).keySet();
			all_nodeIDs_strict.addAll(nodes_strict);
			
			Set<String> nodes_fuzzy = dataAn.getSignificantGenes(data, threshold_fuzzy).keySet();
			nodes_fuzzy.removeAll(nodes_strict);
			all_nodeIDs_fuzzy.addAll(nodes_fuzzy);
			
			System.out.println("  Found " + nodes_strict.size() + " differentially expressed genes at threshold " + threshold_strict + " and " + nodes_fuzzy.size() + " additional ones at threshold " + threshold_fuzzy);
		}
		System.out.println("");
		System.out.println("Defining the set of important nodes");
		
		// Clean out the set of DE genes: if they are strict DE once, they do not need to be in the fuzzy set also
		all_nodeIDs_fuzzy.removeAll(all_nodeIDs_strict);
		System.out.println(" Total: " + all_nodeIDs_strict.size() + " strict differentially expressed genes at threshold " + threshold_strict + " and " + all_nodeIDs_fuzzy.size() + " additional ones at threshold " + threshold_fuzzy);
	
		System.out.println("Expanding the network of DE genes to also include important neighbours");
		
		// Expand the network to include regulatory and PPI partners, and all the connecting fuzzy DE nodes
		Set<String> expandedNetwork = constr.expandNetwork(nm, all_nodeIDs_strict, all_nodeIDs_fuzzy, ppi_file, reg_file, selfInteractions, neighbours, includeUnknownReg);

		// Read all the PPI and regulatory interactions between all the nodes in our expanded network
		// Without modifying edge strenghts, this becomes our reference network
		Set<Node> all_nodes = gp.getNodesByLocusID(expandedNetwork);
		
		System.out.println("Constructing the reference network");
		Set<Edge> ppiEdges = constr.readPPIsByLocustags(nm, ppi_file, all_nodes, all_nodes, selfInteractions);
		System.out.println(" Found " + ppiEdges.size() + " PPI edges between them");
		
		Set<Edge> regEdges = constr.readRegsByLocustags(nm, reg_file, all_nodes, all_nodes, selfInteractions, includeUnknownReg);
		System.out.println(" Found " + regEdges.size() + " PPI regulatory between them");
		
		ppiEdges.addAll(regEdges);
		
		ReferenceNetwork refNet = new ReferenceNetwork("Reference network", firstID++, nm);
		
		// By only defining the edges, unconnected nodes are automatically removed
		refNet.setNodesAndEdges(ppiEdges);
		// refNet.removeUnconnectedNodes();
		
		ReferenceNetwork cleanRefNet = cleaning.fullInputRefCleaning(refNet, nm, eo);
		networks.add(cleanRefNet);
		
		int strictDEnodes1 = 0;
		int fuzzyDEnodes1 = 0;
		
		for (Node n : cleanRefNet.getNodes())
		{
			String id = n.getID();
			if (all_nodeIDs_strict.contains(id))
			{
				strictDEnodes1++;
			}
			if (all_nodeIDs_fuzzy.contains(id))
			{
				fuzzyDEnodes1++;
			}
		}
		System.out.println(" Final, cleaned reference network: " + cleanRefNet.getEdges().size() + " non-redundant edges between " + cleanRefNet.getNodes().size() + 
				" nodes of which " + strictDEnodes1 + " strict DE nodes and " + fuzzyDEnodes1 + " fuzzy DE nodes");

		// Now we create condition-specific networks by altering the edge weights of the original reference network
		for (OverexpressionData data : datasets)
		{
			String suffix = data.getName();
			String name = "Network" + suffix;
			System.out.println("");
			System.out.println("Constructing the condition-specific network for " + name);

			Map<String, Double> all_de_nodes = dataAn.getSignificantGenes(data, threshold_strict);
			all_de_nodes.putAll(dataAn.getSignificantGenes(data, threshold_fuzzy));
			
			Set<Edge> conditionEdges = constr.adjustEdgesByFoldChanges(eo, cleanRefNet.getEdges(), all_de_nodes);
			Condition c = new Condition("time measurement " + suffix);
			ConditionNetwork condNet = new ConditionNetwork(name, firstID++, c, nm);
			
			// By only defining the edges, unconnected nodes are automatically removed
			condNet.setNodesAndEdges(conditionEdges);
			// condNet.removeUnconnectedNodes();
			
			ConditionNetwork cleanCondNet = cleaning.fullInputConditionCleaning(condNet, nm, eo);
			networks.add(cleanCondNet);
			
			int strictDEnodes2 = 0;
			int fuzzyDEnodes2 = 0;
			
			for (Node n : cleanCondNet.getNodes())
			{
				String id = n.getID();
				if (all_nodeIDs_strict.contains(id))
				{
					strictDEnodes2++;
				}
				if (all_nodeIDs_fuzzy.contains(id))
				{
					fuzzyDEnodes2++;
				}
			}
			
			System.out.println(" Condition network " + name + ": " + cleanCondNet.getEdges().size() + " non-redundant edges between " + cleanCondNet.getNodes().size() + 
					" nodes of which " + strictDEnodes2 + " strict DE nodes and " + fuzzyDEnodes2 + " fuzzy DE nodes");
		}
		
		return networks;
	}
	
	/**
	 * Perform the actual 1 against all analysis
	 */
	private void runDiffany(ReferenceNetwork refNet, Set<ConditionNetwork> conditionNets, double weight_cutoff)
	{
		String name = "Osmotic_usecase";
		NodeMapper nm = new DefaultNodeMapper();
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo, nm);
		int runID = p.addRunConfiguration(refNet, conditionNets);
		new CalculateDiff().calculateOneDifferentialNetwork(p, runID, weight_cutoff, 11, 22, true);
		
		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(runID);

		// Testing the edges in the differential network
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
		DifferentialNetwork dNetwork = pair.getDifferentialNetwork();
		
		System.out.println("");
		System.out.println("DifferentialNetwork " + dNetwork + " :");
		for (Edge e : dNetwork.getEdges())
		{
			System.out.println(e);
		}
		
		ConsensusNetwork consNetwork = pair.getConsensusNetwork();
		System.out.println("");
		System.out.println("ConsensusNetwork " + consNetwork + " :");
		for (Edge e : consNetwork.getEdges())
		{
			System.out.println(e);
		}
	}
}
