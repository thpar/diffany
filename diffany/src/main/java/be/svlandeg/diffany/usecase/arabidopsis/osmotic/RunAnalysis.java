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
import be.svlandeg.diffany.core.listeners.ExecutionProgress;
import be.svlandeg.diffany.core.listeners.StandardProgressListener;
import be.svlandeg.diffany.core.networks.Condition;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Node;
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
import be.svlandeg.diffany.usecase.NetworkAnalysis;
import be.svlandeg.diffany.usecase.arabidopsis.ArabidopsisData;
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

	/* Location of all external data. Set to null when you want to exclude either type */
	private URI ppi_file;
	private URI reg_file;
	private URI phos_file;
	private URI kinase_function_file;
	private URI kinase_interaction_file;

	private static String phosAttribute = "phosphorylation_site";
	private static String kinaseAttribute = "kinase_function";

	/**
	 * The constructor defines a few properties of this analysis, such as where to fetch the PPI/regulatory data.
	 */
	public RunAnalysis()
	{
		ppi_file = new ArabidopsisData().getCornetPPI();
		reg_file = new ArabidopsisData().getAtReg();
		phos_file = new ArabidopsisData().getPhosphat();
		kinase_function_file = new ArabidopsisData().getKinases();
		kinase_interaction_file = new ArabidopsisData().getKinaseInteractions();
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

		String inputRoot = "D:" + File.separator + "diffany-osmotic"; // Sofie @ PSB
		//String inputRoot = "C:/Users/Sloffie/Documents/phd/diffany_data/osmotic"; // Sofie @ home

		File osmoticStressDir = new DataIO(inputRoot).getRootOsmoticStressDir();
		String outputDir = osmoticStressDir + File.separator + "output";
		String resultDir = osmoticStressDir + File.separator + "result";

		boolean performStep1FromRaw = false;

		boolean performStep1FromSupplemental = true;
		boolean performStep2ToNetwork = true;
		boolean performStep3InputNetworksToFile = true;

		boolean performStep4InputNetworksFromFile = true;
		boolean performStep5OneagainstAll = false;
		boolean performStep5AllPairwise = true;
		boolean performStep6OutputNetworksToFile = true;

		if (performStep1FromRaw == performStep1FromSupplemental && performStep2ToNetwork)
		{
			System.out.println("Select exactly one option to perform the first step of the analysis!");
			return;
		}
		if (performStep5OneagainstAll && performStep5AllPairwise)
		{
			System.out.println("Select exactly one option to perform the fifth step of the analysis!");
			return;
		}

		if (performStep3InputNetworksToFile && !performStep2ToNetwork)
		{
			System.out.println("Can not perform the third step without the second!");
			return;
		}
		if ((performStep5OneagainstAll || performStep5AllPairwise) && !performStep4InputNetworksFromFile)
		{
			System.out.println("Can not perform the fifth step without the fourth!");
			return;
		}
		if (performStep6OutputNetworksToFile && !(performStep5OneagainstAll || performStep5AllPairwise))
		{
			System.out.println("Can not perform the sixth step without the fifth!");
			return;
		}

		double threshold_strict = 0.05;
		double threshold_fuzzy = 0.1;
		boolean writeHeaders = true;
		boolean allowVirtualEdges = true;

		boolean selfInteractions = false;
		boolean neighbours = true;
		boolean includeUnknownReg = false;
		boolean cleanInputAfterIO = true; // TODO: input should be cleaned before IO in step 3 - but can be set to true to test progresslistener

		double weight_cutoff = 0;
		int hubConnections = 10;
		
		// 5 means that all conditions need to match (4 conditions+reference). 4 means that only 3 time-points (+reference) need to match ("more fuzzy")
		int support = 5;		

		boolean includePhos = true;
		boolean includeKinase = true;
		boolean includePredictedPhos = false;

		String overexpressionFile = null;
		Set<InputNetwork> networks = null;
		RunAnalysis ra = new RunAnalysis();

		/* STEP 1 - OPTION 1: PROCESS RAW CELL DATA TO PRODUCE OVEREXPRESSION VALUES WITH DIFFANY */
		if (performStep1FromRaw)
		{
			System.out.println("1. Transforming CELL data into overexpression values - " + new Date());
			System.out.println("");
			String fileName = "differential_values.txt";
			overexpressionFile = ra.fromRawToOverexpression(osmoticStressDir, fileName);
			System.out.println(" written to " + overexpressionFile);
		}

		/* STEP 1 - OPTION 2: GET PUBLISHED OVEREXPRESSION VALUES FROM THE OSMOTIC PAPER */
		if (performStep1FromSupplemental)
		{
			System.out.println("1. Reading published overexpression values - " + new Date());
			System.out.println("");
			overexpressionFile = osmoticStressDir + File.separator + "clean_Inze_Supplemental_Dataset_1.tab";
		}

		/* STEP 2: USE OVEREXPRESSION VALUES TO CREATE NETWORKS */
		if (performStep2ToNetwork)
		{
			System.out.println("2. Transforming overexpression values into networks - " + new Date());
			System.out.println("   selfinteractions " + selfInteractions + " / neighbours " + neighbours + " / includeUnknownReg " + includeUnknownReg);
			System.out.println("");
			try
			{
				networks = ra.fromOverexpressionToNetworks(new File(overexpressionFile), 1, threshold_strict, threshold_fuzzy, selfInteractions, neighbours,
						includeUnknownReg, includePhos, includeKinase, includePredictedPhos, allowVirtualEdges, hubConnections);
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
		if (performStep2ToNetwork && performStep3InputNetworksToFile)
		{
			System.out.println("");
			System.out.println("3. Writing output networks to " + outputDir + " - " + new Date());
			System.out.println("   allowVirtualEdges " + allowVirtualEdges);
			for (InputNetwork net : networks)
			{
				NetworkIO.writeNetworkToDir(net, net.getNodeMapper(), new File(outputDir, net.getName()), writeHeaders, allowVirtualEdges);
			}
		}

		Set<ConditionNetwork> conditionNets = new HashSet<ConditionNetwork>();
		ReferenceNetwork refNet = null;

		/* STEP 4: READ INPUT NETWORKS BACK IN FROM FILE */
		if (performStep4InputNetworksFromFile)
		{
			System.out.println("");
			System.out.println("4. Reading networks from " + outputDir + " - " + new Date());
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

		/* STEP 5 : GENERATE THE DIFFERENTIAL NETWORKS  */
		RunOutput runOutput = null;
		if (performStep5OneagainstAll || performStep5AllPairwise)
		{
			System.out.println("");
			System.out.println("5. Performing differential analysis at cutoff " + weight_cutoff + " and support " + support + " (only for 1-all) - " + new Date());
			System.out.println("");
			if (refNet == null || conditionNets.isEmpty())
			{
				System.out.println(" Did not find the correct reference and condition-specific networks! ");
				return;
			}
			runOutput = ra.runDiffany(refNet, conditionNets, weight_cutoff, cleanInputAfterIO, support, performStep5AllPairwise);
		}

		/* STEP 6: WRITE THE DIFFERENTIAL AND CONSENSUS NETWORKS TO FILE */
		if (performStep6OutputNetworksToFile)
		{
			System.out.println("");
			System.out.println("6. Writing the generated networks to " + resultDir + " - " + new Date());
			System.out.println("");

			String suffix = "";
			if (performStep5OneagainstAll)
			{
				suffix = "_" + support;
			}

			for (DifferentialNetwork dNetwork : runOutput.getDifferentialNetworks())
			{
				NetworkIO.writeNetworkToDir(dNetwork, dNetwork.getNodeMapper(), new File(resultDir, dNetwork.getName() + suffix), writeHeaders, allowVirtualEdges);
			}

			for (ConsensusNetwork consNetwork : runOutput.getConsensusNetworks())
			{
				NetworkIO.writeNetworkToDir(consNetwork, consNetwork.getNodeMapper(), new File(resultDir, consNetwork.getName() + suffix), writeHeaders, allowVirtualEdges);
			}
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
	private Set<InputNetwork> fromOverexpressionToNetworks(File overExpressionFile, int firstID, double threshold_strict, double threshold_fuzzy,
			boolean selfInteractions, boolean neighbours, boolean includeUnknownReg, boolean includePhos, boolean includeKinase, 
			boolean includePredictedPhos, boolean addVirtualEdges, int hubConnections) throws IOException, URISyntaxException
	{
		Set<String> nodeAttributes = new HashSet<String>();
		if (includePhos)
		{
			nodeAttributes.add(phosAttribute);
		}
		if (includeKinase)
		{
			nodeAttributes.add(kinaseAttribute);
		}

		Set<InputNetwork> networks = new HashSet<InputNetwork>();
		GenePrinter gp = new GenePrinter();

		OverexpressionIO io = new OverexpressionIO();
		ExpressionDataAnalysis dataAn = new ExpressionDataAnalysis();
		NetworkConstruction constr = new NetworkConstruction(gp);

		NodeMapper nm = new DefaultNodeMapper();
		EdgeOntology eo = new DefaultEdgeOntology();
		Logger logger = new Logger();
		NetworkCleaning cleaning = new NetworkCleaning(logger);

		List<OverexpressionData> datasets = io.readDatasets(overExpressionFile, false);

		Set<String> all_nodeIDs_strict = new HashSet<String>();
		Set<String> all_nodeIDs_fuzzy = new HashSet<String>();

		/* Read all different experiments and determine their overexpressed genes */
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

		/* Clean out the set of DE genes: if they are strict DE once, they do not need to be in the fuzzy set also */
		all_nodeIDs_fuzzy.removeAll(all_nodeIDs_strict);
		System.out.println(" Total: " + all_nodeIDs_strict.size() + " strict differentially expressed genes at threshold " + threshold_strict + " and " + all_nodeIDs_fuzzy.size() + " additional ones at threshold " + threshold_fuzzy);

		System.out.println("");
		System.out.println("Expanding the network of DE genes to also include important neighbours");

		/* Expand the network to include regulatory and PPI partners, and all the connecting fuzzy DE nodes */
		Set<String> expandedNetwork = constr.expandNetwork(nm, all_nodeIDs_strict, all_nodeIDs_fuzzy, ppi_file, reg_file, kinase_interaction_file, selfInteractions, neighbours, includeUnknownReg);

		/* Read all the PPI and regulatory interactions between all the nodes in our expanded network
		 * Without modifying edge strenghts, this becomes our reference network */
		Set<Node> all_nodes = gp.getNodesByLocusID(expandedNetwork);

		System.out.println("Constructing the reference network");
		Set<Edge> edges = constr.readPPIsByLocustags(nm, ppi_file, all_nodes, all_nodes, selfInteractions);
		System.out.println(" Found " + edges.size() + " PPI edges between them");

		Set<Edge> regEdges = constr.readRegsByLocustags(nm, reg_file, all_nodes, all_nodes, selfInteractions, includeUnknownReg);
		edges.addAll(regEdges);
		System.out.println(" Found " + regEdges.size() + " regulatory edges between them");
		
		Set<Edge> kinaseEdges = constr.readKinaseInteractionsByLocustags(nm, kinase_interaction_file, all_nodes, all_nodes, selfInteractions);
		edges.addAll(kinaseEdges);
		System.out.println(" Found " + kinaseEdges.size() + " kinase interactions between them");

		// TODO: add as node attributes
		Set<String> phosNodes = constr.readPhosphorylationLocusTags(phos_file, includePredictedPhos);
		Set<String> kinaseNodes = constr.readKinaseLocusTags(kinase_function_file);

		Set<Node> nodes = new HashSet<Node>();
		for (Edge e : edges)
		{
			Node source = e.getSource();
			Node target = e.getTarget();
			nodes.add(source);
			nodes.add(target);
		}
		
		/* Add the appropriate attributes to the nodes */
		for (Node n : nodes)
		{
			if (includePhos)
			{
				if (phosNodes.contains(n.getID()))
				{
					n.setAttribute(phosAttribute, "yes");
				}
				else
				{
					n.setAttribute(phosAttribute, "no");
				}
			}
			if (includeKinase)
			{
				if (kinaseNodes.contains(n.getID()))
				{
					n.setAttribute(kinaseAttribute, "yes");
				}
				else
				{
					n.setAttribute(kinaseAttribute, "no");
				}
			}
		}

		ReferenceNetwork refNet = new ReferenceNetwork("Reference network", firstID++, nodeAttributes, nm);

		/* By only defining the edges, unconnected nodes are automatically removed */
		refNet.setNodesAndEdges(nodes, edges);

		ReferenceNetwork cleanRefNet = cleaning.fullInputRefCleaning(refNet, nm, eo, null);
		networks.add(cleanRefNet);
		nodes = cleanRefNet.getNodes();

		System.out.println(" Final, cleaned reference network: " + cleanRefNet.getEdges().size() + " non-redundant edges between " + cleanRefNet.getNodes().size() + " nodes");

		/* Now we create condition-specific networks by altering the edge weights of the original reference network */
		for (OverexpressionData data : datasets)
		{
			String suffix = data.getName();
			String name = "Network" + suffix;
			Condition c = new Condition("time measurement " + suffix);
			ConditionNetwork condNet = new ConditionNetwork(name, firstID++, nodeAttributes, c, nm);
			
			System.out.println("");
			System.out.println("Constructing the condition-specific network for " + name);

			Map<String, Double> all_de_genes = dataAn.getSignificantGenes(data, threshold_fuzzy);

			NetworkAnalysis na = new NetworkAnalysis();
			String ppiType = "validated_ppi";
			
			Set<String> PPIhubs = na.retrieveHubs(cleanRefNet.getEdges(), nodes, ppiType, hubConnections, false, false);

			Set<Edge> conditionEdges = constr.adjustEdgesByFoldChanges(eo, cleanRefNet.getEdges(), all_de_genes);
			Set<Edge> filteredEdges = constr.filterForHubs(PPIhubs, conditionEdges, ppiType, all_de_genes.keySet());
			
			if (addVirtualEdges)
			{
				Set<Edge> virtualDEedges = constr.constructVirtualRegulations(all_de_genes, nodes);
				
				for (Edge e : virtualDEedges)
				{
					if (includeKinase)
					{
						e.getSource().setAttribute(kinaseAttribute, "unknown");
					}	
					if (includePhos)
					{
						e.getSource().setAttribute(phosAttribute, "unknown");
					}	
				}
				
				filteredEdges.addAll(virtualDEedges);
			}
			
			condNet.setNodesAndEdges(filteredEdges);

			ConditionNetwork cleanCondNet = cleaning.fullInputConditionCleaning(condNet, nm, eo, null);
			networks.add(cleanCondNet);

			System.out.println(" Final, cleaned network " + name + ": " + cleanCondNet.getEdges().size() + " non-redundant edges between " + cleanCondNet.getNodes().size() + " nodes");
		}

		return networks;
	}

	/**
	 * Perform a differential run; either pairwise or 1 against all (no consensus networks)
	 */
	private RunOutput runDiffany(ReferenceNetwork refNet, Set<ConditionNetwork> conditionNets, double weight_cutoff, boolean cleanInputAfterIO, int support, boolean pairwise)
	{
		String name = "Osmotic_usecase_" + support;
		NodeMapper nm = new DefaultNodeMapper();
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo, nm);
		ExecutionProgress listener = new StandardProgressListener();

		int runID = p.addRunConfiguration(refNet, conditionNets, support, cleanInputAfterIO, listener);

		if (pairwise)
		{
			new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, runID, weight_cutoff, true, false, 55, true, listener);
		}
		else
		{
			new CalculateDiff().calculateOneDifferentialNetwork(p, runID, weight_cutoff, 11, -1, true, listener);
		}

		/* Testing that there is exactly one differential network created */
		RunOutput output = p.getOutput(runID);

		return output;
	}
}
