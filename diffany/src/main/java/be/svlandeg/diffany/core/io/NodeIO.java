package be.svlandeg.diffany.core.io;

import java.io.IOException;
import java.util.Set;
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
	 * @return a string representation of all nodes in this network, ready for printing
	 * @see NodeIO#writeToTab
	 */
	public static String writeNodeToTab(Set<Node> nodes)
	{
		String result = "";
		for (Node n : nodes)
		{
			result += NodeIO.writeToTab(n);
			result += System.getProperty("line.separator");
		}
		return result;
	}
	
	/**
	 * Get a string representation of a node.
	 * More specifically, print it as: nodeID - node.name - virtual.
	 * 
	 * @param n the node
	 * @return a string representation of this node, ready for printing
	 */
	public static String writeToTab(Node n)
	{
		String virtualResult = writeVirtualToString(n);
		String result = n.getID() + '\t' + n.getDisplayName() + '\t' + virtualResult;
		return result;
	}
	
	/**
	 * Get a string representation of the boolean virtual state.
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
	 * @return the corresponding Node object
	 * @throws IOException when an error occurs during parsing
	 */
	public static Node readFromTab(String s) throws IOException
	{
		StringTokenizer stok = new StringTokenizer(s, "\t");
		String ID = stok.nextToken();
		String name = stok.nextToken();
		String virtual = stok.nextToken();
		
		boolean isVirtual = NodeIO.isVirtual(virtual);
		
		return new Node(ID, name, isVirtual);
	}
}
