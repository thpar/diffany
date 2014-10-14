package be.svlandeg.diffany.core.networks.meta;

import java.util.Set;

import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.semantics.NodeMapper;

public class MetaConsensusNetwork extends MetaNetwork
{
	
	protected MetaInputNetwork input;
	
	/**
	 * Create a new consensus network, referring to a merged input network.
	 * 
	 * @param name the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network - can be left empty or null
	 * @param nodes the nodes of this network
	 * @param metaEdges the edges of this network
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 * @param input the merged input network
	 * 
	 * @throws IllegalArgumentException when the list of condition-specific networks is null or empty,
	 * or when the reference network is null
	 */
	public MetaConsensusNetwork(String name, int ID, Set<String> nodeAttributes, Set<Node> nodes, Set<MetaEdge> metaEdges, NodeMapper nm, MetaInputNetwork input) 
			throws IllegalArgumentException
	{
		super(name, ID, nodeAttributes, nodes, metaEdges, nm);
		if (input == null)
		{
			String errormsg = "Please define a non-null merged input network!";
			throw new IllegalArgumentException(errormsg);
		}
		this.input = input;
	}
	
	/**
	 * Get the reference network associated to this differential network
	 * @return the reference network 
	 */
	public MetaInputNetwork getInputNetwork()
	{
		return input;
	}
	

	/* (non-Javadoc)
	 * @see be.svlandeg.diffany.core.networks.Network#getStringRepresentation()
	 */
	@Override
	public String getStringRepresentation()
	{
		String result = name + ": merged consensus network calculated from ";
		result += input.getName();
		return result;
	}

}
