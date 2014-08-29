package be.svlandeg.diffany.core.semantics;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.Node;

/**
 * This abstract class allows identifying equal nodes across two networks. 
 * It is not limited to 1-1 node mapping but should be extendable also to N-M mappings.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class NodeMapper
{

	/**
	 * Define whether or not two nodes are equal.
	 * 
	 * @param node1 a node in the first network
	 * @param node2 a node in the second network
	 * @return whether or not they should considered to be equal
	 */
	public abstract boolean areEqual(Node node1, Node node2);

	/**
	 * Return a 'consensus' ID for a set of nodes that were previously determined to be equal. 
	 * Null or empty nodes are ignored. Otherwise, all IDs in the set will be the same, and this will be the result of this method.
	 * 
	 * @param nodes all nodes
	 * @return a consensus name of these nodes to be used in the differential network, or null when all nodes are null
	 * @throws IllegalArgumentException when the nodes are not all equal (when not null)
	 */
	public abstract String getConsensusID(Set<Node> nodes) throws IllegalArgumentException;

	/**
	 * Return a 'consensus' name for a set of nodes that were previously determined to be equal. 
	 * Null or empty nodes are ignored.
	 * 
	 * @param nodes all nodes
	 * @return a consensus name of these nodes to be used in the differential network, or null when all nodes are null
	 * @throws IllegalArgumentException when the nodes are not all equal (when not null)
	 */
	public abstract String getConsensusName(Set<Node> nodes) throws IllegalArgumentException;

	/**
	 * Define whether or not this node is already in the set.
	 * 
	 * @param node a node to be tested for its presence in the set
	 * @param nodeSet a set of non-redundant nodes
	 * @return whether or not the set already includes (an equal of) the node
	 */
	public boolean isContained(Node node, Set<Node> nodeSet)
	{
		boolean contained = false;
		for (Node n : nodeSet)
		{
			if (areEqual(n, node))
			{
				contained = true;
			}
		}
		return contained;
	}

	/**
	 * Define all equal nodes in the two networks, mapping one node in network 1 to a set of nodes in network 2.
	 * 
	 * @param network1 the first network
	 * @param network2 the second network
	 * @return all equal nodes, mapping network 1 to network 2
	 */
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

	/**
	 * Define all equal nodes in a set of networks, creating disjoint sets of nodes.
	 * 
	 * @param networks all input networks
	 * @return all equal nodes grouped in sets
	 */
	public List<Set<Node>> getAllEquals(Set<Network> networks)
	{
		List<Set<Node>> allEqualSets = new ArrayList<Set<Node>>();

		for (Network n : networks)
		{
			// for each network, try to find the equals of each node
			for (Node queryNode : n.getNodes())
			{
				boolean foundEquals = false;
				for (Set<Node> equalsSet : allEqualSets)
				{
					Node compareNode = equalsSet.iterator().next();
					if (!foundEquals && compareNode != null)
					{
						// we assume equality transitivity: if one node is equal to our querynode, then all will be
						if (areEqual(queryNode, compareNode))
						{
							foundEquals = true;
							equalsSet.add(queryNode);
						}
					}
				}
				if (! foundEquals)
				{
					Set<Node> newEqualsSet = new HashSet<Node>();
					newEqualsSet.add(queryNode);
					allEqualSets.add(newEqualsSet);
				}
			}
		}
		return allEqualSets;
	}

	/**
	 * Get all nodes, without duplicating equal nodes. 
	 * 
	 * @param networks the networks
	 * @return the union of all nodes in all networks, removing duplicates using the areEqual method
	 */
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

}
