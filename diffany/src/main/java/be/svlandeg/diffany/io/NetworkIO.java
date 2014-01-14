package be.svlandeg.diffany.io;

import java.util.Set;

import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.Network;

/**
 * This class allows reading or writing a {@link Network} from File.
 * 
 * @author Sofie Van Landeghem
 */
public class NetworkIO
{
	
	/**
	 * Get a string representation of all edges in a Network, divided by newlines, with edges in a tabbed format.
	 * 
	 * @return a string representation of all edges in this network, ready for printing
	 * @see EdgeIO#writeToTab
	 */
	public static String writeEdgesToTab(Network n)
	{
		String result = "";
		for (Edge e : n.getEdges())
		{
			result += EdgeIO.writeToTab(e);
			result += System.getProperty("line.separator");
		}
		return result;
	}
	
	/**
	 * Read all edges from a file containing one tab-delimited edge per line
	 * 
	 * @param location the file location of the tab-delimited file
	 * @return the set of edges read from the file
	 */
	public static Set<Edge> readEdgesFromFile(String location)
	{
		//TODO v1.1: implement!
		throw new UnsupportedOperationException("readEdgesFromFile not yet implemented");
	}
	
	

}
