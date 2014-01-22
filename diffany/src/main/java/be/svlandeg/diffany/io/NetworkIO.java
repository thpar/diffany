package be.svlandeg.diffany.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
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
	 * @param n the Network that needs to be written
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
	 * Write the string representation of all edges to a File.
	 * 
	 * @param f the output file
	 * @param n the Network that needs to be written
	 * @see EdgeIO#writeToTab
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static void writeNetworkToFile(Network n, File f) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(f));
		for (Edge e : n.getEdges())
		{
			writer.append(EdgeIO.writeToTab(e));
			writer.newLine();
			writer.flush();
		}
		writer.flush();
		writer.close();
	}

	/**
	 * Read all edges from a file containing one tab-delimited edge per line
	 * 
	 * @param location the file location of the tab-delimited file
	 * @return the set of edges read from the file
	 */
	public static Set<Edge> readEdgesFromFile(String location)
	{
		// TODO v1.1: implement!
		throw new UnsupportedOperationException("readEdgesFromFile not yet implemented");
	}

}
