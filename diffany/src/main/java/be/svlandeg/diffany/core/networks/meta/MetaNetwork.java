package be.svlandeg.diffany.core.networks.meta;

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
public abstract class MetaNetwork extends Network
{
	/**
	 * Create a new network with a specific name and sets of nodes and edges.
	 * All source and target nodes of each edge will be automatically added to the internal set of nodes.
	 * 
	 * @param name the name of this network 
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network - can be left empty or null
	 * @param nodes the nodes of this network
	 * @param metaEdges the edges of this network
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 */
	public MetaNetwork(String name, int ID, Set<String> nodeAttributes, Set<Node> nodes, Set<MetaEdge> metaEdges, NodeMapper nm)
	{
		super(name, ID, nodeAttributes, nm);
		if (nm == null)
		{
			String errormsg = "Please define a proper NodeMapper object!";
			throw new IllegalArgumentException(errormsg);
		}

		this.name = name;
		this.nm = nm;
		Set<Edge> edges = MetaConvertor.castToNormalEdges(metaEdges);
		setNodesAndEdges(nodes, edges);
	}

	/**
	 * Create a new network with an empty set of nodes and edges.
	 * 
	 * @param name the name of this network 
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network - can be left empty or null
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 */
	public MetaNetwork(String name, int ID, Set<String> nodeAttributes, NodeMapper nm)
	{
		this(name, ID, nodeAttributes, new HashSet<Node>(), new HashSet<MetaEdge>(), nm);
	}

}
