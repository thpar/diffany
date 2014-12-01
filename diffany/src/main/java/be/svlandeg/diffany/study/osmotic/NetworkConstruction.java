package be.svlandeg.diffany.study.osmotic;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.EdgeDefinition;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.semantics.EdgeOntology;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.study.osmotic.arabidopsis.KinaseData;
import be.svlandeg.diffany.study.osmotic.arabidopsis.PPIdata;
import be.svlandeg.diffany.study.osmotic.arabidopsis.RegData;

/**
 * This class allows to construct networks out of overexpression/coexpression values.
 * 
 * @author Sofie Van Landeghem
 */
public class NetworkConstruction
{

	// TODO v2.1 make this class more generic / generalizable

	private GenePrinter gp;

	/**
	 * Constructor: defines the gene printer that can deal with gene synonymy etc.
	 * 
	 * @param gp the gene printer object
	 */
	public NetworkConstruction(GenePrinter gp)
	{
		this.gp = gp;
	}

	/**
	 * This method defines all the nodes that will be in the Diffany networks to analyse a given set of overexpressed genes.
	 * The set of strict and 'fuzzy' DE genes thus contain all genes that are DE in at least one of the conditions in the experiment.
	 * 
	 * @param nodeIDs_strict_DE the overexpressed nodes (not null), with stringent criteria (e.g. FDR 0.05)
	 * @param nodeIDs_fuzzy_DE the overexpressed nodes (may be null), with less stringent criteria (e.g. FDR 0.1, this set does not need to include the strict DE ones)
	 * @param selfInteractions whether or not to include self interactions
	 * @param neighbours whether or not to include direct neighbours
	 * @param includeUnknownReg whether or not to include unknown regulations
	 * 
	 * @return the set of found edges
	 * @throws IOException when an IO error occurs
	 * @throws URISyntaxException when an input file could not be read
	 */
	public Set<String> expandNetwork(Set<String> nodeIDs_strict_DE, Set<String> nodeIDs_fuzzy_DE, boolean selfInteractions, boolean neighbours, boolean includeUnknownReg) throws URISyntaxException, IOException
	{
		NetworkAnalysis na = new NetworkAnalysis();
		// TODO v2.1: this method unnecessarily generates Network objects?

		Set<String> allNodeIDs = new HashSet<String>();
		allNodeIDs.addAll(nodeIDs_strict_DE);
		Set<Node> nodes_strict_DE = gp.getNodesByLocusID(nodeIDs_strict_DE);

		PPIdata ppi = new PPIdata(gp);
		KinaseData kinase = new KinaseData(gp);
		RegData reg = new RegData(gp);

		/* 1. expand the (strict) DE node set with PPI neighbours */

		/* First define the PPI hubs */
		Set<Edge> allAthPPIedges = ppi.readAllPPIs(selfInteractions);
		Network ppiAthNetwork_all = new InputNetwork("Test network PPI", 341, null);
		ppiAthNetwork_all.setNodesAndEdges(allAthPPIedges);
		Set<String> PPIhubs = na.retrieveHubs(allAthPPIedges, ppiAthNetwork_all.getNodes(), PPIdata.hubPPI, false, true, true);
		Set<Node> nodes_PPI_hubs = gp.getNodesByLocusID(PPIhubs);

		Set<Edge> PPIedges_strict = null;
		if (neighbours)
		{
			/* Here, we exclude PPI hubs as non-DE neighbours and also exclude them as DE partners that extend the network */
			PPIedges_strict = ppi.readPPIsByLocustags(nodes_strict_DE, nodes_PPI_hubs, null, nodes_PPI_hubs, selfInteractions);
		}
		else
		{
			PPIedges_strict = ppi.readPPIsByLocustags(nodes_strict_DE, null, nodes_strict_DE, null, selfInteractions);
		}

		InputNetwork PPInetwork_strict = new InputNetwork("PPI network", 342, null);
		PPInetwork_strict.setNodesAndEdges(PPIedges_strict);
		allNodeIDs.addAll(NodeMapper.getNodeIDs(PPInetwork_strict.getNodes()));

		/* 2. expand the original (strict) DE node set with regulatory and kinase neighbours */
		Set<Edge> regEdges1 = new HashSet<Edge>();
		if (neighbours)
		{
			regEdges1.addAll(reg.readRegsByLocustags(nodes_strict_DE, null, selfInteractions, includeUnknownReg)); // from our input to their targets
			regEdges1.addAll(reg.readRegsByLocustags(null, nodes_strict_DE, selfInteractions, includeUnknownReg)); // from our input to their sources (may result in redundant nodes but these will be cleaned out later)
		}

		if (neighbours)
		{
			/* First define the hubs */
			Set<Edge> allAthKinaseEdges = kinase.readAllKinaseInteractions(selfInteractions);
			Network kinaseAthNetwork = new InputNetwork("Test network kinase", 343, null);
			kinaseAthNetwork.setNodesAndEdges(allAthKinaseEdges);
			Set<String> kinaseTargetHubs = na.retrieveHubs(allAthKinaseEdges, kinaseAthNetwork.getNodes(), KinaseData.hubPhos, false, true, false);
			Set<String> kinaseSourceHubs = na.retrieveHubs(allAthKinaseEdges, kinaseAthNetwork.getNodes(), KinaseData.hubPhos, false, false, true);
			Set<Node> nodes_kinase_source_hubs = gp.getNodesByLocusID(kinaseSourceHubs);
			Set<Node> nodes_kinase_target_hubs = gp.getNodesByLocusID(kinaseTargetHubs);

			regEdges1.addAll(kinase.readKinaseInteractionsByLocustags(nodes_strict_DE, nodes_kinase_source_hubs, null, nodes_kinase_target_hubs, selfInteractions)); // from our input to their targets
			regEdges1.addAll(kinase.readKinaseInteractionsByLocustags(null, nodes_kinase_source_hubs, nodes_strict_DE, nodes_kinase_target_hubs, selfInteractions)); // from our input to their sources (may result in redundant nodes but these will be cleaned out later)
		}

		InputNetwork regNetwork1 = new InputNetwork("Regulatory network", 666, null);
		regNetwork1.setNodesAndEdges(regEdges1);
		Set<String> expandedNodeIDs2 = NodeMapper.getNodeIDs(regNetwork1.getNodes());

		allNodeIDs.addAll(expandedNodeIDs2);

		/* 3. add all fuzzy DE nodes which connect to the strict DE nodes or the PPI/regulatory partners */
		if (nodeIDs_fuzzy_DE != null && !nodeIDs_fuzzy_DE.isEmpty())
		{
			Set<Node> allNodes = gp.getNodesByLocusID(allNodeIDs);
			Set<Node> nodes_fuzzy_DE = gp.getNodesByLocusID(nodeIDs_fuzzy_DE);

			Set<Edge> PPIedges_fuzzy = ppi.readPPIsByLocustags(allNodes, null, nodes_fuzzy_DE, null, selfInteractions);
			InputNetwork PPInetwork = new InputNetwork("PPI network", 342, null);
			PPInetwork.setNodesAndEdges(PPIedges_fuzzy);
			Set<String> expandedPPINodeIDs = NodeMapper.getNodeIDs(PPInetwork.getNodes());

			Set<Edge> regEdges2 = new HashSet<Edge>();
			
			regEdges2.addAll(reg.readRegsByLocustags(nodes_fuzzy_DE, allNodes, selfInteractions, includeUnknownReg)); // from fuzzy DE to our combined set
			regEdges2.addAll(reg.readRegsByLocustags(allNodes, nodes_fuzzy_DE, selfInteractions, includeUnknownReg)); // from our combined set to fuzzy DE (may result in redundant nodes but these will be cleaned out later)

			regEdges2.addAll(kinase.readKinaseInteractionsByLocustags(nodes_fuzzy_DE, null, allNodes, null, selfInteractions)); // from fuzzy DE to our combined set
			regEdges2.addAll(kinase.readKinaseInteractionsByLocustags(allNodes, null, nodes_fuzzy_DE, null, selfInteractions)); // from our combined set to fuzzy DE (may result in redundant nodes but these will be cleaned out later)

			InputNetwork regNetwork2 = new InputNetwork("Regulatory network", 666, null);
			regNetwork2.setNodesAndEdges(regEdges2);
			Set<String> expandedRegNodeIDs = NodeMapper.getNodeIDs(regNetwork2.getNodes());

			allNodeIDs.addAll(expandedPPINodeIDs);
			allNodeIDs.addAll(expandedRegNodeIDs);
		}

		return allNodeIDs;
	}

	/**
	 * Take a copy of the given set of edges and produce new ones that are weighted according to the fold changes of the DE genes.
	 * 
	 * @param eo the edge ontology which can tell which edge types are anti-correlated (e.g. inhibition)
	 * @param nodes the set of nodes that should be used to construct the new edges
	 * @param origEdges the original set of edges (will not be changed!)
	 * @param all_de_nodes all overexpressed nodes (both stringent and less-stringent criteria)
	 * @return a new set of edges, whose weights are changed according to the fold changes of the DE genes
	 */
	public Set<Edge> adjustEdgesByFoldChanges(EdgeOntology eo, Set<Node> nodes, Set<Edge> origEdges, Map<String, Double> all_de_nodes)
	{
		Map<String, Node> mappedNodes = NodeMapper.getNodesByID(nodes);

		Set<Edge> resultEdges = new HashSet<Edge>();
		for (Edge e : origEdges)
		{
			String sourceID = e.getSource().getID();
			String targetID = e.getTarget().getID();

			double sourceFC = 0;
			if (all_de_nodes.containsKey(sourceID))
			{
				sourceFC = all_de_nodes.get(sourceID);
			}

			double targetFC = 0;
			if (all_de_nodes.containsKey(targetID))
			{
				targetFC = all_de_nodes.get(targetID);
			}

			boolean correlation = true;
			if (!e.isSymmetrical())
			{
				String edgeType = e.getType();
				String edgeCat = eo.getSourceCategory(edgeType);
				if (eo.getAllNegSourceCategories().contains(edgeCat))
				{
					// this means they are anti-correlated, e.g. by negative regulation
					correlation = false;
				}
			}

			double weight = 1;

			// both are non-DE (we need to keep this if-statement for the following else statements to work properly)
			if (targetFC == 0 && sourceFC == 0)
			{
				weight = 1;
			}

			// both are up-regulated or non-DE: the positive edge weight is (1 + their FC average), or 0 in case they are anti-correlated
			else if (targetFC >= 0 && sourceFC >= 0)
			{
				weight = 1 + ((sourceFC + targetFC) / 2);
				if (!correlation && targetFC != 0) // if target is non-DE, the weight can remain positive
				{
					weight = 0;
				}
			}

			// the target is downregulated (below 0), the source is up (above 0) or non-DE -> break the connection unless they are anti-correlated
			else if (targetFC < 0 && sourceFC >= 0)
			{
				weight = 0;
				if (!correlation)
				{
					weight = 1 + ((sourceFC + Math.abs(targetFC)) / 2);
				}
			}

			// the source is downregulated, the target is upregulated or non-DE -> break the connection, even if they are anti-correlated
			// if they were anti-correlated (e.g. inhibition), this is not the case anymore and the differential will be an increase_regulation
			else if (sourceFC < 0 && targetFC >= 0)
			{
				weight = 0;
			}

			// both are down-regulated: the positive edge weight ]0,1[ is their FC average through an exponential function for normalization, 
			// or 0 if they are anti-correlated
			else if (targetFC < 0 && sourceFC < 0)
			{
				weight = 0;
				if (correlation)
				{
					double neg_avg = (sourceFC + targetFC) / 2;

					// this 0.5 factor determines the steepness of the normalization curve
					// after this normalization, the weight will be somewhere between 0 and 1
					weight = Math.exp(0.5 * neg_avg);
				}
			}

			EdgeDefinition newEdge = new EdgeDefinition(e.getDefinition());
			newEdge.setWeight(weight);
			resultEdges.add(new Edge(mappedNodes.get(sourceID), mappedNodes.get(targetID), newEdge));
		}

		return resultEdges;
	}

	/**
	 * Modify the differentially expressed status of a collection of DE nodes, using over/under expression data.
	 * 
	 * @param DEnodeIDs the map of target node IDs with their corresponding fold change values (negative weights refer to down-regulation)
	 * @param nodes the set of existing nodes in the network - they will acquire additional attributes after this method
	 */
	public void modifyDEState(Map<String, Double> DEnodeIDs, Set<Node> nodes)
	{
		String de_attribute = Node.de_attribute;
		for (Node n : nodes)
		{
			String id = n.getID();
			Double FC = DEnodeIDs.get(id);
			if (FC == null || FC == 0)
			{
				n.setAttribute(de_attribute, Node.not_de);
			}
			else if (FC < 0)
			{
				n.setAttribute(de_attribute, Node.downregulated);
			}
			else if (FC > 0)
			{
				n.setAttribute(de_attribute, Node.upregulated);
			}
		}
	}

}
