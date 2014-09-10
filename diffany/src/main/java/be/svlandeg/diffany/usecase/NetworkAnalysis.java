package be.svlandeg.diffany.usecase;

import java.net.URI;
import java.util.HashMap;
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
	 * @param n the input network
	 */
	public void analyseHubs(Network n)
	{
		System.out.println("HUB ANALYSIS");
		Map<String, Integer> edgeCountByNodeID = new HashMap<String, Integer>();
		for (Node node : n.getNodes())
		{
			String id = node.getID();
			edgeCountByNodeID.put(id, 0);
		}
		
		for (Edge e : n.getEdges())
		{
			String sourceID = e.getSource().getID();
			Integer count1 = edgeCountByNodeID.get(sourceID);
			edgeCountByNodeID.put(sourceID, count1+1);
			
			String targetID = e.getSource().getID();
			Integer count2 = edgeCountByNodeID.get(targetID);
			edgeCountByNodeID.put(targetID, count2+1);
		}
		
		SortedMap<Integer, Integer> edgeCountByOccurrence = new TreeMap<Integer, Integer>();
		for (Node node : n.getNodes())
		{
			String id = node.getID();
			Integer edgeCount = edgeCountByNodeID.get(id);
			Integer occurrence = edgeCountByOccurrence.get(edgeCount);
			if (occurrence == null)
			{
				occurrence = 0;
			}
			edgeCountByOccurrence.put(edgeCount, occurrence+1);
			
			
		}
		int totalOccEdges = 0;
		int totalOccNodes = 0;
		for (Integer edgeCount : edgeCountByOccurrence.keySet())
		{
			int occurrence = edgeCountByOccurrence.get(edgeCount);
			totalOccEdges += edgeCount * occurrence;
			totalOccNodes += occurrence;
		}
		int accumulOccEdges = 0;
		int accumulOccNodes = 0;
		for (Integer edgeCount : edgeCountByOccurrence.keySet())
		{
			int occurrence = edgeCountByOccurrence.get(edgeCount);
			System.out.println(" - " + edgeCount + " edges -> " + occurrence + " occurrences (unique nodes)");
			accumulOccEdges += edgeCount * occurrence;
			accumulOccNodes += occurrence;
			double percEdges = 100 * accumulOccEdges / totalOccEdges;
			double percNodes = 100 * accumulOccNodes / totalOccNodes;
			System.out.println("   Accumulated occurrence: " + accumulOccEdges + "/" + totalOccEdges + " = " + percEdges + "% edges for " 
					+ accumulOccNodes + "/" + totalOccNodes + " = " + percNodes + "% nodes");
		}
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
		na.analyseHubs(network);
		
		System.out.println(" ");
		System.out.println("DONE");
	}

}
