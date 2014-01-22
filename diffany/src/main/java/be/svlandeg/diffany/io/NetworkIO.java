package be.svlandeg.diffany.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Node;

/**
 * This class allows reading or writing a {@link Network} from File.
 * 
 * @author Sofie Van Landeghem
 */
public class NetworkIO
{

	private static String default_edge_file = "edges.tab";
	private static String default_node_file = "nodes.tab";
	
	
	/**
	 * Write the string representation of the network: All edges to one File (edges.tab), and all nodes to another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files will be written
	 * @param network the Network that needs to be written
	 * @see EdgeIO#writeToTab
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	public static void writeNetworkToDir(Network network, File dir) throws IOException
	{
		File edgeFile = new File(dir.getAbsolutePath() + "/" + default_edge_file);
		File nodeFile = new File(dir.getAbsolutePath() + "/" + default_node_file);
		writeNetworkToFiles(network, edgeFile, nodeFile);
	}

	/**
	 * Write the string representation of the network: All edges to one File, and all nodes to another, as specified by the parameters
	 * 
	 * @param edgesFile the output file in which the edges will be written
	 * @param nodesFile the output file in which the edges will be written
	 * @param network the Network that needs to be written
	 * @see EdgeIO#writeToTab
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	public static void writeNetworkToFiles(Network network, File edgesFile, File nodesFile) throws IOException
	{
		BufferedWriter edgeWriter = new BufferedWriter(new FileWriter(edgesFile));
		for (Edge e : network.getEdges())
		{
			edgeWriter.append(EdgeIO.writeToTab(e));
			edgeWriter.newLine();
			edgeWriter.flush();
		}
		edgeWriter.flush();
		edgeWriter.close();

		BufferedWriter nodeWriter = new BufferedWriter(new FileWriter(nodesFile));
		for (Node n : network.getNodes())
		{
			nodeWriter.append(n.getName());
			nodeWriter.newLine();
			nodeWriter.flush();
		}
		nodeWriter.flush();
		nodeWriter.close();
	}
	
	/**
	 * Read a network from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param name the name of the network (TODO: should this be written to file as well?)
	 * @return a Network representation of the nodes and edges in the files.
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static Network readNetworkFromDir(File dir, String name) throws IOException
	{
		File edgeFile = new File(dir.getAbsolutePath() + "/" + default_edge_file);
		File nodeFile = new File(dir.getAbsolutePath() + "/" + default_node_file);
		
		Set<Edge> edges = readEdgesFromFile(edgeFile);
		Set<Node> nodes = readNodesFromFile(nodeFile);
		
		// TODO: decide on the type of network? Or allow an abstract network?
		return null;
	}

	/**
	 * Read all edges from a file containing one tab-delimited edge per line
	 * 
	 * @param edgesFile the file containing the edge data
	 * @return the set of edges read from the file
	 * @see EdgeIO#readFromTab
	 * 
	 * @throws IOException when an errors occurs during parsing
	 */
	public static Set<Edge> readEdgesFromFile(File edgesFile) throws IOException
	{
		Set<Edge> edges = new HashSet<Edge>();

		BufferedReader reader = new BufferedReader(new FileReader(edgesFile));
		String line = reader.readLine();
		while (line != null)
		{
			Edge e = EdgeIO.readFromTab(line);
			edges.add(e);
			line = reader.readLine();
		}

		reader.close();

		return edges;
	}
	
	/**
	 * Read all nodes from a file containing one node name per line
	 * 
	 * @param nodesFile the file containing the node data
	 * @return the set of nodes read from the file
	 * 
	 * @throws IOException when an errors occurs during parsing
	 */
	public static Set<Node> readNodesFromFile(File nodesFile) throws IOException
	{
		Set<Node> nodes = new HashSet<Node>();

		BufferedReader reader = new BufferedReader(new FileReader(nodesFile));
		String line = reader.readLine();
		
		while (line != null)
		{
			// TODO : should we somehow check the line is not ridiculously long?
			Node n = new Node(line);
			nodes.add(n);
			line = reader.readLine();
		}

		reader.close();
		return nodes;
	}

}
