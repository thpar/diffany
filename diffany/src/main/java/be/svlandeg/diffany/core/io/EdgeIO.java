package be.svlandeg.diffany.core.io;

import java.io.IOException;
import java.util.Set;
import java.util.StringTokenizer;

import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.EdgeDefinition;
import be.svlandeg.diffany.core.networks.Node;

/**
 * This class allows reading/writing an {@link Edge} from/to File.
 * 
 * @author Sofie Van Landeghem
 */
public class EdgeIO
{
	
	private static String symmString = "symmetrical";
	private static String directString = "directed";
	private static String negatedString = "negated";
	private static String notnegatedString = "not negated";

	/**
	 * Get a string representation of all edges in a collection, divided by newlines, with edges in a tabbed format.
	 * 
	 * @param edges the edges that needs to be written
	 * @return a string representation of all edges in this network, ready for printing
	 * @see EdgeIO#writeToTab
	 */
	public static String writeEdgesToTab(Set<Edge> edges)
	{
		String result = "";
		for (Edge e : edges)
		{
			result += EdgeIO.writeToTab(e);
			result += System.getProperty("line.separator");
		}
		return result;
	}
	
	/**
	 * Get a string representation of an edge.
	 * More specifically, print it as: source.name - target.name - edge.type - symmetrical - weight - negated.
	 * 
	 * @param e the edge
	 * @return a string representation of this edge, ready for printing
	 */
	public static String writeToTab(Edge e)
	{
		String defResult = writeDefinitionToTab(e);
		String result = e.getSource().getName() + '\t' + e.getTarget().getName() + '\t' + defResult;
		return result;
	}
	
	
	/**
	 * Read an Edge from a tab-delimited String.
	 * 
	 * @param s the tab-delimited string containing all parameters of the Edge
	 * @return the corresponding Edge object
	 * @throws IOException when an error occurs during parsing
	 */
	public static Edge readFromTab(String s) throws IOException
	{
		StringTokenizer stok = new StringTokenizer(s, "\t");
		String source = stok.nextToken();
		String target = stok.nextToken();
		
		String defS = "";
		while (stok.hasMoreTokens())
		{
			defS += stok.nextToken() + "\t";
		}
		
		EdgeDefinition def =  readDefinitionToTab(defS);
		return new Edge(new Node(source), new Node(target), def);
	}
	

	/**
	 * Get a string representation of an edge definition.
	 * More specifically, print it as: edge.type - symmetrical - weight - negated.
	 * 
	 * @return a string representation of this edge definition, ready for printing
	 */
	public static String writeDefinitionToTab(EdgeDefinition def)
	{
		String symm = symmString;
		if (! def.isSymmetrical())
		{
			symm = directString;
		}
		String neg = negatedString;
		if (! def.isNegated())
		{
			neg = notnegatedString;
		}
		double weight = def.getWeight();
		String result = def.getType() + '\t' + symm + '\t' + weight + '\t' + neg;
		
		return result;
	}
	
	
	/**
	 * Read an EdgeDefinition from a tab-delimited String.
	 * 
	 * @return a string representation of this edge , ready for printing
	 * @throws IOException when an error occurs during parsing
	 */
	public static EdgeDefinition readDefinitionToTab(String def) throws IOException
	{
		StringTokenizer stok = new StringTokenizer(def, "\t");
		String type = stok.nextToken();
		String symm_string = stok.nextToken();
		double weight = Double.parseDouble(stok.nextToken());
		String neg_string = stok.nextToken();
		
		boolean symmetrical = true;
		if (symm_string.equals(directString))
		{
			symmetrical = false;
		}
		else if (! symm_string.equals(symmString))
		{
			throw new IOException("Error reading the symmetry state from the edge :" + symm_string);
		}
		
		boolean negated = true;
		if (neg_string.equals(notnegatedString))
		{
			negated = false;
		}
		else if (! neg_string.equals(negatedString))
		{
			throw new IOException("Error reading the negation state from the edge :" + neg_string);
		}
		return new EdgeDefinition(type, symmetrical, weight, negated);
	}
	
	

}
