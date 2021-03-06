package be.svlandeg.diffany.study.osmotic;

/*
 * #%L
 * Diffany
 * %%
 * Copyright (C) 2014 PSB/UGent - Sofie Van Landeghem and Thomas Van Parys
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU General Lesser Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/lgpl-3.0.html>.
 * #L%
 */


import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.algorithms.NetworkCleaning;
import be.svlandeg.diffany.core.io.NetworkIO;
import be.svlandeg.diffany.core.networks.Attribute;
import be.svlandeg.diffany.core.networks.Condition;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.progress.ProgressListener;
import be.svlandeg.diffany.core.progress.StandardProgressListener;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.EdgeOntology;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;
import be.svlandeg.diffany.study.osmotic.arabidopsis.GenePrinter;
import be.svlandeg.diffany.study.osmotic.arabidopsis.KinaseData;
import be.svlandeg.diffany.study.osmotic.arabidopsis.PPIdata;
import be.svlandeg.diffany.study.osmotic.arabidopsis.RegData;

/**
 * This class provides the analysis pipeline for the osmotic stress application on Arabidopsis. 
 * The main method calls all procedures necessary to perform the full analysis, from data input to result output.
 * 
 * @author Sofie Van Landeghem
 */
public class RunAnalysis
{

	/**
	 * Empty constructor.
	 */
	public RunAnalysis()
	{}

	/**
	 * Run the full analysis pipeline.
	 * The input directory should be given as argument. 
	 * This is where the files 'input' and 'result' will be written/read, and where a file 'clean_Inze_Supplemental_Dataset_1.tab' should be located.
	 * 
	 * @param args CL arguments - the first should contain the input directory
	 * @throws Exception whenever something goes wrong
	 */
	public static void main(String[] args) throws Exception
	{
		if (args.length < 1)
		{
			System.out.println("The first argument should point to the input directory!");
			return;
		}

		String inputLocation = args[0] + File.separator; 
		
		System.out.println("Performing osmotic data analysis from " + inputLocation + " - " + new Date());
		System.out.println("");
		
		String outputDir = inputLocation + "output";
		String resultDir = inputLocation + "result";
		
		boolean performStep1FromSupplemental = false;
		boolean performStep2ToNetwork = false;
		boolean performStep3InputNetworksToFile = false;

		boolean performStep4InputNetworksFromFile = true;
		boolean performStep5OneagainstAll = true;
		boolean performStep5AllPairwise = false;
		boolean performStep6OutputNetworksToFile = false;

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
		double threshold_fuzzy = 0.1;		// if you don't want a fuzzy cut-off, simply put it equal to the strict one.
		
		if (threshold_fuzzy < threshold_strict)
		{
			System.out.println("threshold_fuzzy " + threshold_fuzzy + " needs to be higher than threshold_strict " + threshold_strict + "!");
			return;
		}
		
		boolean writeHeaders = true;

		boolean selfInteractions = false;
		boolean neighbours = true;
		boolean includeUnknownReg = false;
		
		// the input should have been cleaned before IO in step 3 - but this can be set to true to test e.g. the progresslistener or cleaning speed
		boolean cleanInputAfterIO = false; 

		double weight_cutoff = 0;
		
		// e.g. 5 means that all conditions need to match (4 conditions+reference). 4 means that only 3 time-points (+reference) need to match ("more fuzzy")
		int support = 5;		

		boolean includePhos = true;
		boolean includeKinase = true;
		boolean includePredictedPhos = false;

		String overexpressionFile = null;
		Set<InputNetwork> networks = null;
		RunAnalysis ra = new RunAnalysis();


		/* STEP 1: GET PUBLISHED OVEREXPRESSION VALUES FROM THE OSMOTIC PAPER */
		if (performStep1FromSupplemental)
		{
			System.out.println("1. Reading published overexpression values - " + new Date());
			System.out.println("");
			overexpressionFile = inputLocation + "clean_Inze_Supplemental_Dataset_1.tab";
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
						includeUnknownReg, includePhos, includeKinase, includePredictedPhos);
			}
			catch (IllegalArgumentException e)
			{
				System.out.println("Trouble parsing from " + overexpressionFile + ": " + e.getMessage());
				e.printStackTrace();
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
			for (InputNetwork net : networks)
			{
				NetworkIO.writeNetworkToDir(net, new File(outputDir, net.getName()), writeHeaders);
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

			Set<InputNetwork> readNetworks = NetworkIO.readGenericInputNetworksFromSubdirs(new File(outputDir), writeHeaders, true);
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
				NetworkIO.writeNetworkToDir(dNetwork, new File(resultDir, dNetwork.getName() + suffix), writeHeaders);
			}

			for (ConsensusNetwork consNetwork : runOutput.getConsensusNetworks())
			{
				NetworkIO.writeNetworkToDir(consNetwork, new File(resultDir, consNetwork.getName() + suffix), writeHeaders);
			}
		}

		System.out.println("");
		System.out.println("Done! - " + new Date());
	}


	/**
	 * Use the overexpression values to generate input networks: 1 reference network + 1 condition-dependent network for each overexpression dataset that can be read from overExpressionFile
	 */
	private Set<InputNetwork> fromOverexpressionToNetworks(File overExpressionFile, int firstID, double threshold_strict, double threshold_fuzzy,
			boolean selfInteractions, boolean neighbours, boolean includeUnknownReg, boolean includePhos, boolean includeKinase, 
			boolean includePredictedPhos) throws IOException, URISyntaxException
	{
		Set<Attribute> nodeAttributes = new HashSet<Attribute>();
		nodeAttributes.add(new Attribute(Node.de_attribute, String.class));
		if (includePhos)
		{
			nodeAttributes.add(new Attribute(Node.phos_attribute, String.class));
		}
		if (includeKinase)
		{
			nodeAttributes.add(new Attribute(Node.kinase_attribute, String.class));
		}

		Set<InputNetwork> networks = new HashSet<InputNetwork>();
		GenePrinter gp = new GenePrinter();

		OverexpressionIO io = new OverexpressionIO();
		NetworkConstruction constr = new NetworkConstruction(gp);
		
		PPIdata ppi = new PPIdata(gp);
		KinaseData kinase = new KinaseData(gp);
		RegData reg = new RegData(gp);

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

			Set<String> ids_strict = data.getSignificantGenes(threshold_strict).keySet();
			all_nodeIDs_strict.addAll(ids_strict);

			Set<String> ids_fuzzy = data.getSignificantGenes(threshold_fuzzy).keySet();
			ids_fuzzy.removeAll(ids_strict);
			all_nodeIDs_fuzzy.addAll(ids_fuzzy);

			System.out.println("  Found " + ids_strict.size() + " differentially expressed genes at threshold " + threshold_strict + " and " + ids_fuzzy.size() + " additional ones at threshold " + threshold_fuzzy);
		}
		System.out.println("");
		System.out.println("Defining the set of important nodes");

		/* Clean out the set of DE genes: if they are strict DE once, they do not need to be in the fuzzy set also */
		all_nodeIDs_fuzzy.removeAll(all_nodeIDs_strict);
		System.out.println(" Total: " + all_nodeIDs_strict.size() + " strict differentially expressed genes at threshold " + threshold_strict + " and " + all_nodeIDs_fuzzy.size() + " additional ones at threshold " + threshold_fuzzy);

		System.out.println("");
		System.out.println("Expanding the network of DE genes to also include important neighbours");

		/* Expand the network to include regulatory and PPI partners, and all the connecting fuzzy DE nodes */
		Set<String> expanded_ID_set = constr.expandNetwork(all_nodeIDs_strict, all_nodeIDs_fuzzy, selfInteractions, neighbours, includeUnknownReg);

		/* Read all the PPI and regulatory interactions between all the nodes in our expanded network
		 * Without modifying edge strenghts, this becomes our reference network */
		Set<Node> ref_nodes = gp.getNodesByLocusID(expanded_ID_set);
		
		System.out.println("Constructing the reference network");
		Set<Edge> ref_edges = ppi.readPPIsByLocustags(ref_nodes, null, ref_nodes, null, selfInteractions);
		System.out.println(" Found " + ref_edges.size() + " PPI edges between them");

		Set<Edge> ref_reg_Edges = reg.readRegsByLocustags(ref_nodes, ref_nodes, selfInteractions, includeUnknownReg);
		System.out.println(" Found " + ref_reg_Edges.size() + " regulatory edges between them");
		ref_edges.addAll(ref_reg_Edges);
		
		Set<Edge> ref_kinase_edges = kinase.readKinaseInteractionsByLocustags(ref_nodes, null, ref_nodes, null, selfInteractions);
		System.out.println(" Found " + ref_kinase_edges.size() + " kinase interactions between them");
		ref_edges.addAll(ref_kinase_edges);

		Set<String> phosNodes = kinase.readPhosphorylationLocusTags(includePredictedPhos);
		Set<String> kinaseNodes = kinase.readKinaseLocusTags();

		Set<Node> ref_filtered_nodes = new HashSet<Node>();
		for (Edge e : ref_edges)
		{
			Node source = e.getSource();
			Node target = e.getTarget();
			ref_filtered_nodes.add(source);
			ref_filtered_nodes.add(target);
		}
		
		/* Add the appropriate attributes to the nodes */
		for (Node n : ref_filtered_nodes)
		{
			setFixedAttributes(n, phosNodes, kinaseNodes, includePhos, includeKinase);
			n.setAttribute(Node.de_attribute, Node.not_de);
		}

		ReferenceNetwork refNet = new ReferenceNetwork("Reference network", firstID++, nodeAttributes);

		/* By only defining the edges, unconnected nodes are automatically removed */
		refNet.setNodesAndEdges(ref_filtered_nodes, ref_edges);

		ReferenceNetwork cleanRefNet = cleaning.fullInputRefCleaning(refNet, eo, null);
		networks.add(cleanRefNet);
		
		Set<String> ref_node_IDs = new HashSet<String>();
		for (Node n : cleanRefNet.getNodes())
		{
			ref_node_IDs.add(n.getID());
		}

		System.out.println(" Final, cleaned reference network: " + cleanRefNet.getEdges().size() + " non-redundant edges between " + cleanRefNet.getNodes().size() + " nodes");

		/* Now we create condition-specific networks by altering the edge weights of the original reference network */
		for (OverexpressionData data : datasets)
		{
			// create a new copy of nodes
			Set<Node> condition_nodes = gp.getNodesByLocusID(ref_node_IDs);
			
			String suffix = data.getName();
			String name = "Network" + suffix;
			Condition c = new Condition("time measurement " + suffix);
			ConditionNetwork condNet = new ConditionNetwork(name, firstID++, nodeAttributes, c);
			
			System.out.println("");
			System.out.println("Constructing the condition-specific network for " + name);

			// record the DE state as node attribute
			Map<String, Double> all_de_genes = data.getSignificantGenes(threshold_fuzzy);
			constr.modifyDEState(all_de_genes, condition_nodes);
			
			for (Node n : condition_nodes)
			{
				setFixedAttributes(n, phosNodes, kinaseNodes, includePhos, includeKinase);
			}

			Set<Edge> condition_edges = constr.adjustEdgesByFoldChanges(eo, condition_nodes, cleanRefNet.getEdges(), all_de_genes);
			condNet.setNodesAndEdges(condition_edges);
			ConditionNetwork cleanCondNet = cleaning.fullInputConditionCleaning(condNet, eo, null);
			networks.add(cleanCondNet);

			System.out.println(" Final, cleaned network " + name + ": " + cleanCondNet.getEdges().size() + " non-redundant edges between " + cleanCondNet.getNodes().size() + " nodes");
		}

		return networks;
	}
	
	/**
	 * Set the node attributes for phosphorylation site and kinase function.
	 */
	private void setFixedAttributes(Node n, Set<String> phosNodes, Set<String> kinaseNodes, boolean includePhos, boolean includeKinase)
	{
		if (includePhos)
		{
			if (phosNodes.contains(n.getID()))
			{
				n.setAttribute(Node.phos_attribute, "yes");
			}
			else
			{
				n.setAttribute(Node.phos_attribute, "no");
			}
		}
		if (includeKinase)
		{
			if (kinaseNodes.contains(n.getID()))
			{
				n.setAttribute(Node.kinase_attribute, "yes");
			}
			else
			{
				n.setAttribute(Node.kinase_attribute, "no");
			}
		}
	}

	/**
	 * Perform a differential run; either pairwise or 1 against all (no consensus networks)
	 */
	private RunOutput runDiffany(ReferenceNetwork refNet, Set<ConditionNetwork> conditionNets, double weight_cutoff, boolean cleanInputAfterIO, int support, boolean pairwise)
	{
		String name = "Osmotic_usecase_" + support;
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo);
		ProgressListener listener = new StandardProgressListener(true);

		int runID = p.addRunConfiguration(refNet, conditionNets, support, cleanInputAfterIO, listener);

		if (pairwise)
		{
			new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, runID, weight_cutoff, true, false, 55, true, listener);
		}
		else
		{
			new CalculateDiff().calculateOneDifferentialNetwork(p, runID, weight_cutoff, null, null, 11, -1, true, listener);
		}
		RunOutput output = p.getOutput(runID);
		return output;
	}
}
