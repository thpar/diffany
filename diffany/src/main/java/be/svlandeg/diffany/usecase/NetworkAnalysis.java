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
	 * Print information about hubs in the network. This method disregards symmetry/directionality information, but simply counts all edges, once for source, and one for the target.
	 * This method is always applied for one specific type of interactions, e.g. 'ppi', in order not to mix different types of connections.
	 * Further, this method returns a set of edges which has min_connections or more.
	 * 
	 * @param edges the edges in the input network 
	 * @param nodes the nodes in the input network 
	 * @param type the type of interactions we want to count hubs for
	 * @param min_connections the number of connections a 'hub' needs to have to be defined as such (inclusive)
	 * @param printDetails whether or not to print a detailed account of the hub analysis
	 * @param divideByTwo put this parameter to true if symmetrical edges are represented twice in the data (i.e. the edges are not cleaned)
	 * @return the set of node IDs which are considered to be hubs at the specified percentage cutoff
	 */
	public Set<String> retrieveHubs(Set<Edge> edges, Set<Node> nodes, String type, int min_connections, boolean printDetails, boolean divideByTwo)
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
			// only count this edge if it is of the correct type
			if (e.getType().equals(type))
			{
				String sourceID = e.getSource().getID();
				Integer count1 = edgeCountByNodeID.get(sourceID);
				edgeCountByNodeID.put(sourceID, count1+1);
				
				String targetID = e.getSource().getID();
				Integer count2 = edgeCountByNodeID.get(targetID);
				edgeCountByNodeID.put(targetID, count2+1);
			}
		}
		
		SortedMap<Integer, Set<String>> occurrenceByEdgeCount = new TreeMap<Integer, Set<String>>();
		for (Node node : nodes)
		{
			String id = node.getID();
			Integer edgeCount = edgeCountByNodeID.get(id);
			if (divideByTwo)
			{
				edgeCount = edgeCount / 2;
			}
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
		int connections = 10;
		Set<String> PPIhubs = na.retrieveHubs(network.getEdges(), network.getNodes(), "validated_ppi", connections, printDetails, true);
		System.out.println(" Found " + PPIhubs.size() + " PPI hubs:");
		for (String hub : PPIhubs)
		{
			System.out.println("  " + hub);
		}
		
		System.out.println(" ");
		System.out.println("DONE");
	}

}
