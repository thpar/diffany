package be.svlandeg.diffany.core.networks.merged;

import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.semantics.NodeMapper;

/**
 * A kind of network that merges information from different conditions, using edges which know which conditions they belong to.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class MergedNetwork extends Network
{
	/**
	 * Create a new network with a specific name and sets of nodes and edges.
	 * All source and target nodes of each edge will be automatically added to the internal set of nodes.
	 * 
	 * @param name the name of this network (should be enforced to be unique within one project)
	 * @param nodes the nodes of this network
	 * @param mergedEdges the edges of this network
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 */
	public MergedNetwork(String name, Set<Node> nodes, Set<MergedEdge> mergedEdges, NodeMapper nm)
	{
		super(name, nm);
		if (nm == null)
		{
			String errormsg = "Please define a proper NodeMapper object!";
			throw new IllegalArgumentException(errormsg);
		}

		this.name = name;
		this.nm = nm;
		Set<Edge> edges = MergedConvertor.castToNormalEdges(mergedEdges);
		setNodesAndEdges(nodes, edges);
	}

	/**
	 * Create a new network with an empty set of nodes and edges.
	 * @param name the name of this network (should be enforced to be unique within one project)
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 */
	public MergedNetwork(String name, NodeMapper nm)
	{
		this(name, new HashSet<Node>(), new HashSet<MergedEdge>(), nm);
	}

}
