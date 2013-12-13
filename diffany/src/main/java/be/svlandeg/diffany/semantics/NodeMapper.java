package be.svlandeg.diffany.semantics;

import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Node;

/**
 * This interface allows identifying equal nodes across two networks. 
 * It's implementation may be project-specific, for instance a CyNodeMapper may expect specifically CyNodes (or such), 
 * as the NodeMapper is defined at a time when the networks are known (when making a new Project).
 * 
 * @author Sofie Van Landeghem
 */
public interface NodeMapper
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
	 * Return a 'consensus' name for a set of nodes that were previously determined to be equal. 
	 * Null or empty nodes are ignored.
	 * 
	 * @param nodes all nodes
	 * @return a consensus name of these nodes to be used in the differential network, or null when all nodes are null
	 * @throws IllegalArgumentException when the nodes are not all equal (when not null)
	 */
	public abstract String getConsensusName(Set<Node> nodes) throws IllegalArgumentException;

	/**
	 * Define all equal nodes in the two networks, mapping one node in network 1 to a set of nodes in network 2.
	 * 
	 * @param network1 the first network
	 * @param network2 the second network
	 * @return all equal nodes, mapping network 1 to network 2
	 */
	public abstract Map<Node, Set<Node>> getAllEquals(Network network1, Network network2);
	
	/**
	 * Get all nodes, without duplicating equal nodes. 
	 * 
	 * @param networks the networks
	 * @return the union of all nodes in all networks, removing duplicates using the areEqual method
	 */
	public abstract Set<Node> getAllNodes(Set<Network> networks);

}
