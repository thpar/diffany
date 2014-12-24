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
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.StringTokenizer;

import be.svlandeg.diffany.core.networks.Attribute;
import be.svlandeg.diffany.core.networks.Node;

/**
 * This class allows reading/writing a {@link Node} from/to File.
 * 
 * @author Sofie Van Landeghem
 */
public class NodeIO
{

	
	/**
	 * Get a string representation of all nodes in a collection, divided by newlines, with nodes in a tabbed format.
	 * 
	 * @param nodes the nodes that needs to be written
	 * @param nodeAttributes the node attribute names
	 * 
	 * @return a string representation of all nodes in this network, ready for printing
	 * @see NodeIO#writeToTab
	 */
	public static String writeNodesToTab(Set<Node> nodes, SortedSet<Attribute> nodeAttributes)
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
	 * More specifically, print it as: nodeID - node.name - and then the node attributes (alphabetically).
	 * 
	 * @param n the node
	 * @param nodeAttributes the node attribute names
	 * 
	 * @return a string representation of this node, ready for printing
	 */
	public static String writeToTab(Node n, SortedSet<Attribute> nodeAttributes)
	{
		String result = n.getID() + '\t' + n.getDisplayName();
		for (Attribute att : nodeAttributes)
		{
			Object value = n.getAttribute(att.getName());
			String stringValue = new String();
			if (value != null){
				stringValue = value.toString();
			}
			result += '\t' + stringValue;
		}
		return result;
	}
	
	/**
	 * Write the header for the node file, including the custom node attribute names.
	 * @param nodeAttributes the node attributes defined for the network of which this node is a part of
	 * @return a String containing the tab-delimited header for the node file
	 */
	public static CharSequence getHeader(SortedSet<Attribute> nodeAttributes)
    {
		String result = "ID" + '\t' + "official_symbol";
		for (Attribute attribute : nodeAttributes)
		{
			result += "\t" + attribute;
		}
		return result;
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
		
		Node n = new Node(ID, name);
		
		for (String att : nodeAttributes)
		{
			String value = stok.nextToken();
			n.setAttribute(att, value);
		}
		
		return n;
	}

}
