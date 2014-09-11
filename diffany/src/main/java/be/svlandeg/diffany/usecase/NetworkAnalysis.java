package be.svlandeg.diffany.usecase;

import java.net.URI;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.usecase.arabidopsis.CornetData;
import be.svlandeg.diffany.usecase.arabidopsis.GenePrinter;
import be.svlandeg.diffany.usecase.arabidopsis.NetworkConstruction;


/**
 * This class is useful for generating general statistics about Diffany networks.
 * 
 * @author Sofie Van Landeghem
 */
public class NetworkAnalysis
{
	
	/**
	 * Currently empty constructor
	 */
	public NetworkAnalysis()
	{
		
	}
	
	/**
	 * Print a few generic network stats
	 * @param n the input network
	 */
	public void printNetworkStats(Network n)
	{
		System.out.println("NETWORK ANALYSIS");
		System.out.println(" - " + n.getEdges().size() + " edges");
		System.out.println(" - " + n.getNodes().size() + " nodes");
	}
	
	/**
	 * Print information about hubs in the network. This method currently disregards symmetry/directionality information (TODO).
	 * 
	 * @param edges the edges in the input network 
	 * @param nodes the nodes in the input network 
	 * @param perc_most_connected the percentage of most connected hubs we want to obtain
	 * @param printDetails whether or not to print a detailed account of the hub analysis
	 * @return the set of node IDs which are considered to be hubs at the specified percentage cutoff
	 */
	public Set<String> analyseHubs(Set<Edge> edges, Set<Node> nodes, int perc_most_connected, boolean printDetails)
	{
 		Set<String> hubs = new HashSet<String>();
		Map<String, Integer> edgeCountByNodeID = new HashMap<String, Integer>();
		for (Node node : nodes)
		{
			String id = node.getID();
			edgeCountByNodeID.put(id, 0);
		}
		
		for (Edge e : edges)
		{
			String sourceID = e.getSource().getID();
			Integer count1 = edgeCountByNodeID.get(sourceID);
			edgeCountByNodeID.put(sourceID, count1+1);
			
			String targetID = e.getSource().getID();
			Integer count2 = edgeCountByNodeID.get(targetID);
			edgeCountByNodeID.put(targetID, count2+1);
		}
		
		SortedMap<Integer, Set<String>> occurrenceByEdgeCount = new TreeMap<Integer, Set<String>>();
		for (Node node : nodes)
		{
			String id = node.getID();
			Integer edgeCount = edgeCountByNodeID.get(id);
			Set<String> occurrences = occurrenceByEdgeCount.get(edgeCount);
			if (occurrences == null)
			{
				occurrences = new HashSet<String>();
			}
			occurrences.add(id);
			occurrenceByEdgeCount.put(edgeCount, occurrences);
			
			
		}
		int totalOccEdges = 0;
		int totalOccNodes = 0;
		for (Integer edgeCount : occurrenceByEdgeCount.keySet())
		{
			int occurrence = occurrenceByEdgeCount.get(edgeCount).size();
			totalOccEdges += edgeCount * occurrence;
			totalOccNodes += occurrence;
		}
		int accumulOccEdges = 0;
		int accumulOccNodes = 0;
		for (Integer edgeCount : occurrenceByEdgeCount.keySet())
		{
			Set<String> occurrences = occurrenceByEdgeCount.get(edgeCount);
			
			accumulOccEdges += edgeCount * occurrences.size();
			accumulOccNodes += occurrences.size();
			double percEdges = 100 * accumulOccEdges / totalOccEdges;
			double percNodes = 100 * accumulOccNodes / totalOccNodes;
			
			if (percNodes >= perc_most_connected)
			{
				hubs.addAll(occurrences);
			}
			
			if (printDetails)
			{
				System.out.println(" - " + edgeCount + " edges -> " + occurrences.size() + " occurrences (unique nodes)");
				System.out.println("   Accumulated occurrence: " + accumulOccEdges + "/" + totalOccEdges + " = " + percEdges + "% edges for " 
					+ accumulOccNodes + "/" + totalOccNodes + " = " + percNodes + "% nodes");
			}
		}
		return hubs;
	}
	
	
	/**
	 * Run an analysis on specific data
	 * @param args run parameters (which are ignored)
	 * @throws Exception when an IO error or something such occurs
	 */
	public static void main(String[] args) throws Exception
	{
		NetworkAnalysis na = new NetworkAnalysis();
		GenePrinter gp = new GenePrinter();
		NetworkConstruction constr = new NetworkConstruction(gp);
		NodeMapper nm = new DefaultNodeMapper();
		
		boolean includeSelfInteractions = false;
		URI file = new CornetData().getCornetPPI();
		
		System.out.println("Reading: " + file + " - " + "includeSelfInteractions=" + includeSelfInteractions);
		
		Set<Edge> edges = constr.readAllPPIs(nm, file, includeSelfInteractions);
		Network network = new InputNetwork("Test network", 342, nm);
		network.setNodesAndEdges(edges);
		
		System.out.println(" ");
		na.printNetworkStats(network);
		
		System.out.println(" ");
		System.out.println("HUB ANALYSIS");
		boolean printDetails = true;
		int perc = 98;
		Set<String> hubs = na.analyseHubs(network.getEdges(), network.getNodes(), perc, printDetails);
		System.out.println(" Found " + hubs.size() + " hubs:");
		for (String hub : hubs)
		{
			System.out.println("  " + hub);
		}
		
		System.out.println(" ");
		System.out.println("DONE");
	}

}
