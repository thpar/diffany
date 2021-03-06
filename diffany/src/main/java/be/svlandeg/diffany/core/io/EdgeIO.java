package be.svlandeg.diffany.core.io;

/*
 * #%L
 * Diffany
 * %%
 * Copyright (C) 2014 PSB/UGent - Sofie Van Landeghem and Thomas Van Parys
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU General Lesser Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/lgpl-3.0.html>.
 * #L%
 */


import java.io.IOException;
import java.util.Map;
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
	private static String directString = "asymmetrical";
	
	private static String negatedString = "negated";
	private static String notnegatedString = "affirmative";

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
			result += writeToTab(e);
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
		String defResult = writeDefinitionToTab(e.getDefinition());
		String result = e.getSource().getID() + '\t' + e.getTarget().getID() + '\t' + defResult;
		return result;
	}
	
	/**
	 * Get a string representation of an edge definition.
	 * More specifically, print it as: edge.type - symmetrical - weight - negated.
	 * 
	 * @param def the original edge definition
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
	 * Write the header line, i.e. a tab-delimited summary of the information that is printed with writeToTab.
	 * @return a tab-delimited string respresentation of edge data, useful as header in text files
	 */
	public static String getHeader()
    {
		String result = "source" + '\t' + "target" + '\t' + "interaction" + '\t' + "symmetry" + '\t' + "weight" + '\t' + "negation";
		return result;
    }
	
	
	/**
	 * Read an Edge from a tab-delimited String. A set of input nodes should be given as input and mapped by their unique ID,
	 * to be able to link the edge to the correct nodes.
	 * 
	 * @param s the tab-delimited string containing all parameters of the Edge
	 * @param nodes the input nodes in the network, containing at least the two nodes for this edge (otherwise source/target will be null)
	 * @return the corresponding Edge object
	 * @throws IOException when an error occurs during parsing
	 */
	public static Edge readFromTab(String s, Map<String, Node> nodes) throws IOException
	{
		StringTokenizer stok = new StringTokenizer(s, "\t");
		String source = stok.nextToken();
		String target = stok.nextToken();
		
		String defS = "";
		while (stok.hasMoreTokens())
		{
			defS += stok.nextToken() + "\t";
		}
		
		EdgeDefinition def =  readDefinitionFromTab(defS);
		return new Edge(nodes.get(source), nodes.get(target), def);
	}
	
	
	
	/**
	 * Read an EdgeDefinition from a tab-delimited String.
	 * 
	 * @param def the original edge definition, in string format
	 * @return the edge definition represented by the input string
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static EdgeDefinition readDefinitionFromTab(String def) throws IOException
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
