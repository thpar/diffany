package be.svlandeg.diffany.core.networks.merged;

import java.util.Set;

import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.semantics.NodeMapper;

/**
 * A kind of network that merges information from different conditions, using edges which know which conditions they belong to.
 * This specific class is used for the input data.
 * 
 * @author Sofie Van Landeghem
 */
public class MergedInputNetwork extends MergedNetwork
{
	
	/**
	 * Create a new input network with a specific name and sets of nodes and edges.
	 * All source and target nodes of each edge will be automatically added to the internal set of nodes.
	 * 
	 * @param name the name of this network (should be enforced to be unique within one project)
	 * @param nodes the nodes of this network
	 * @param conditionEdges the edges of this network
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 */
	public MergedInputNetwork(String name, Set<Node> nodes, Set<ConditionEdge> conditionEdges, NodeMapper nm)
	{
		super(name, nodes, conditionEdges, nm);
	}

	@Override
	public String getStringRepresentation()
	{
		return name + ": merged input network";
	}

}
