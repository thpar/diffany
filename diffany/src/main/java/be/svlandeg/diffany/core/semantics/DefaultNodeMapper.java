package be.svlandeg.diffany.core.semantics;

import java.util.*;

import be.svlandeg.diffany.core.networks.Node;

/**
 * This class provides a default implementation of a NodeMapper, defining nodes as equal when their IDs are equal.
 * 
 * @author Sofie Van Landeghem
 */
public class DefaultNodeMapper extends NodeMapper
{

	@Override
	public boolean areEqual(Node node1, Node node2)
	{
		if (node1 == null && node2 != null)
			return false;
		if (node1 == null && node2 == null)
			return true;
		return node1.getID().equals(node2.getID());
	}

	
	@Override
	public String getConsensusID(Set<Node> nodes) throws IllegalArgumentException
	{
		Set<String> IDs = new HashSet<String>();
		for (Node n : nodes)
		{
			if (n != null)
			{
				IDs.add(n.getID());
			}
		}
		if (IDs.isEmpty())
		{
			return null;
		}
		if (IDs.size() > 1)
		{
			String found = "";
			for (String s : IDs)
			{
				found += s + " / ";
			}
			String errormsg = "A consensus ID can only be defined for two equal nodes (same IDs), but found more than one ID: " + found;
			throw new IllegalArgumentException(errormsg);
		}
		
		return IDs.iterator().next();
	}
}
