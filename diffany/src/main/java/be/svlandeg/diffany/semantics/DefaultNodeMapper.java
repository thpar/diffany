package be.svlandeg.diffany.semantics;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

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
	public String getConsensusName(Node node1, Node node2) throws IllegalArgumentException
	{
		if (!areEqual(node1, node2))
		{
			String errormsg = "A consensus name can only be defined for two equal nodes!";
			throw new IllegalArgumentException(errormsg);
		}
		return node1.getName(true);
	}
}
