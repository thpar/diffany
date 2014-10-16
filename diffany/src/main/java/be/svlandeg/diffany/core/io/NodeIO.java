package be.svlandeg.diffany.core.io;

import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.StringTokenizer;

import be.svlandeg.diffany.core.networks.Node;

/**
 * This class allows reading/writing a {@link Node} from/to File.
 * 
 * @author Sofie Van Landeghem
 */
public class NodeIO
{
	
	private static String virtualString = "virtual";
	private static String nonVirtualString = "normal";

	
	/**
	 * Get a string representation of all nodes in a collection, divided by newlines, with nodes in a tabbed format.
	 * 
	 * @param nodes the nodes that needs to be written
	 * @param nodeAttributes the node attribute names
	 * 
	 * @return a string representation of all nodes in this network, ready for printing
	 * @see NodeIO#writeToTab
	 */
	public static String writeNodesToTab(Set<Node> nodes, SortedSet<String> nodeAttributes)
	{
		String result = "";
		for (Node n : nodes)
		{
			result += NodeIO.writeToTab(n, nodeAttributes);
			result += System.getProperty("line.separator");
		}
		return result;
	}
	
	/**
	 * Get a string representation of a node.
	 * More specifically, print it as: nodeID - node.name - virtual and then the node attributes (alphabetically).
	 * 
	 * @param n the node
	 * @param nodeAttributes the node attribute names
	 * 
	 * @return a string representation of this node, ready for printing
	 */
	public static String writeToTab(Node n, SortedSet<String> nodeAttributes)
	{
		String virtualResult = writeVirtualToString(n);
		String result = n.getID() + '\t' + n.getDisplayName() + '\t' + virtualResult;
		for (String att : nodeAttributes)
		{
			result += '\t' + n.getAttribute(att);
		}
		return result;
	}
	
	/**
	 * Write the header for the node file, including the custom node attribute names.
	 * @param nodeAttributes the node attributes defined for the network of which this node is a part of
	 * @return a String containing the tab-delimited header for the node file
	 */
	public static CharSequence getHeader(SortedSet<String> nodeAttributes)
    {
		String result = "ID" + '\t' + "official_symbol" + '\t' + "virtual";
		for (String attribute : nodeAttributes)
		{
			result += "\t" + attribute;
		}
		return result;
    }
	
	/**
	 * Get a string representation of the boolean virtual state.
	 * @param n the original node
	 * 
	 * @return a string representation of this node's virtual state
	 */
	public static String writeVirtualToString(Node n)
	{
		String result = nonVirtualString;
		if (n.isVirtual())
		{
			result = virtualString;
		}
		return result;
	}
	
	
	/**
	 * Parse the boolean virtual state from the Node's string representation.
	 * 
	 * @param s the string representation
	 * @return a boolean indicating the virtual state of the node
	 */
	public static boolean isVirtual(String s)
	{
		if (s.equals(virtualString))
		{
			return true;
		}
		return false;
	}
	
	/**
	 * Read a Node from a tab-delimited String.
	 * 
	 * @param s the tab-delimited string containing all parameters of the Node
	 * @param nodeAttributes the node attribute names
	 * @return the corresponding Node object
	 * @throws IOException when an error occurs during parsing
	 */
	public static Node readFromTab(String s, List<String> nodeAttributes) throws IOException
	{
		StringTokenizer stok = new StringTokenizer(s, "\t");
		String ID = stok.nextToken();
		String name = stok.nextToken();
		String virtual = stok.nextToken();
		boolean isVirtual = NodeIO.isVirtual(virtual);
		
		Node n = new Node(ID, name, isVirtual);
		
		for (String att : nodeAttributes)
		{
			String value = stok.nextToken();
			n.setAttribute(att, value);
		}
		
		return n;
	}

}
