package be.svlandeg.diffany.semantics;

import java.util.*;

import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Node;

/**
 * This class provides a default implementation of a NodeMapper, defining nodes
 * as equal when there names are equal (ignoring case, i.e. using 'normalized'
 * names).
 * 
 * @author Sofie Van Landeghem
 * 
 */
public class DefaultNodeMapper implements NodeMapper
{

	@Override
	public boolean areEqual(Node node1, Node node2)
	{
		if (node1 == null && node2 != null)
			return false;
		if (node1 == null && node2 == null)
			return true;
		return node1.getName(true).equals(node2.getName(true));
	}

	@Override
	public Map<Node, Set<Node>> getAllEquals(Network network1, Network network2)
	{
		Map<Node, Set<Node>> allEquals = new HashMap<Node, Set<Node>>();
		for (Node node1 : network1.getNodes())
		{
			Set<Node> equalNodes = new HashSet<Node>();
			for (Node node2 : network2.getNodes())
			{
				if (areEqual(node1, node2))
				{
					equalNodes.add(node2);
				}
			}
			allEquals.put(node1, equalNodes);
		}
		return allEquals;
	}

	@Override
	public Set<Node> getAllNodes(Set<Network> networks)
	{
		Set<Node> allNodes = new HashSet<Node>();
		for (Network network : networks)
		{
			for (Node node : network.getNodes())
			{
				boolean hasEqual = false;
				for (Node compNode : allNodes)
				{
					if (areEqual(node, compNode))
					{
						hasEqual = true;
					}
				}
				if (!hasEqual)
				{
					allNodes.add(node);
				}
			}
		}
		return allNodes;
	}

	@Override
	public String getConsensusName(Set<Node> nodes) throws IllegalArgumentException
	{
		Set<String> normalized_names = new HashSet<String>();
		Set<String> original_names = new HashSet<String>();
		for (Node n : nodes)
		{
			if (n != null)
			{
				original_names.add(n.getName(false));
				normalized_names.add(n.getName(true));
			}
		}
		if (original_names.isEmpty())
		{
			return null;
		}
		if (normalized_names.size() > 1)
		{
			String found = "";
			for (String s : normalized_names)
			{
				found += s + " / ";
			}
			String errormsg = "A consensus name can only be defined for two equal nodes, found: " + found;
			throw new IllegalArgumentException(errormsg);
		}
		if (original_names.size() == 1)
		{
			// keep casing information if it is the same
			return original_names.iterator().next();
		}
		// return lowercase otherwise
		return normalized_names.iterator().next();
	}

}
