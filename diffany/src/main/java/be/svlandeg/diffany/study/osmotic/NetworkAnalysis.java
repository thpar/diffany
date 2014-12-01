package be.svlandeg.diffany.study.osmotic;

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
	 * Print information about hubs in the network. This method disregards symmetry/directionality information, but simply counts all edges, once for source, and one for the target.
	 * This method is typically applied for one specific type of interactions, but in principle the given set of edges can also be mixed.
	 * Further, this method returns a set of edges which has min_connections or more.
	 * 
	 * @param edges the edges in the input network 
	 * @param nodes the nodes in the input network 
	 * @param min_connections the number of connections a 'hub' needs to have to be defined as such (inclusive)
	 * @param printDetails whether or not to print a detailed account of the hub analysis
	 * @param countIncoming whether or not to count outgoing edges (when dealing with symmetrical edges, incoming and outgoing should both be true)
	 * @param countOutgoing whether or not to count incoming edges (when dealing with symmetrical edges, incoming and outgoing should both be true)
	 * @return the set of node IDs which are considered to be hubs at the specified percentage cutoff
	 */
	public Set<String> retrieveHubs(Set<Edge> edges, Set<Node> nodes, int min_connections, boolean printDetails, boolean countIncoming, boolean countOutgoing)
	{
		if (printDetails)
		{
			System.out.println("edge count \t unique nodes \t accumulated edges \t percEdges \t accumulated nodes \t percNodes");
		}
		
 		Set<String> hubs = new HashSet<String>();
		Map<String, Integer> edgeCountByNodeID = new HashMap<String, Integer>();
		for (Node node : nodes)
		{
			String id = node.getID();
			edgeCountByNodeID.put(id, 0);
		}
		
		for (Edge e : edges)
		{
			if (countOutgoing)
			{
				String sourceID = e.getSource().getID();
				Integer count1 = edgeCountByNodeID.get(sourceID);
				edgeCountByNodeID.put(sourceID, count1+1);
			}
			
			if (countIncoming)
			{
				String targetID = e.getTarget().getID();
				Integer count2 = edgeCountByNodeID.get(targetID);
				edgeCountByNodeID.put(targetID, count2+1);
			}
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
		for (int edgeCount = 0; edgeCount <= occurrenceByEdgeCount.lastKey(); edgeCount++)
		{
			Set<String> occurrences = occurrenceByEdgeCount.get(edgeCount);
			if (occurrences == null)
			{
				occurrences = new HashSet<String>();
			}
			
			accumulOccEdges += edgeCount * occurrences.size();
			accumulOccNodes += occurrences.size();
			
			int percEdges = 100 * accumulOccEdges / totalOccEdges;
			
			int percNodes = 100 * accumulOccNodes / totalOccNodes;
			
			if (edgeCount >= min_connections)
			{
				hubs.addAll(occurrences);
			}
			
			if (printDetails)
			{
				System.out.println(edgeCount + "\t" + occurrences.size() + "\t" + accumulOccEdges + "\t" + percEdges + "\t" + accumulOccNodes + "\t" + percNodes);
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
		
		boolean includeSelfInteractions = false;
		URI file = new ArabidopsisData().getCornetPPI();
		
		System.out.println("Reading: " + file + " - " + "includeSelfInteractions=" + includeSelfInteractions);
		
		Set<Edge> PPIedges = constr.readAllPPIs(file, includeSelfInteractions);
		Network network = new InputNetwork("Test network", 342, null);
		network.setNodesAndEdges(PPIedges);
		
		System.out.println(" ");
		na.printNetworkStats(network);
		
		System.out.println(" ");
		System.out.println("HUB ANALYSIS");
		boolean printDetails = true;
		int connections = 10;
		Set<String> PPIhubs = na.retrieveHubs(PPIedges, network.getNodes(), connections, printDetails, true, true);
		System.out.println(" Found " + PPIhubs.size() + " PPI hubs:");
		for (String hub : PPIhubs)
		{
			System.out.println("  " + hub);
		}
		
		System.out.println(" ");
		System.out.println("DONE");
	}

}
