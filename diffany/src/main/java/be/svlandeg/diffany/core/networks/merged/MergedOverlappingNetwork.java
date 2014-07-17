package be.svlandeg.diffany.core.networks.merged;

import java.util.Set;

import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.semantics.NodeMapper;

public class MergedOverlappingNetwork extends MergedNetwork
{
	
	protected MergedInputNetwork input;
	
	/**
	 * Create a new overlapping network, referring to a merged input network.
	 * 
	 * @param name the name of this network
	 * @param nodes the nodes of this network
	 * @param mergedEdges the edges of this network
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 * @param input the merged input network
	 * 
	 * @throws IllegalArgumentException when the list of condition-specific networks is null or empty,
	 * or when the reference network is null
	 */
	public MergedOverlappingNetwork(String name, Set<Node> nodes, Set<MergedEdge> mergedEdges, NodeMapper nm, MergedInputNetwork input) 
			throws IllegalArgumentException
	{
		super(name, nm);
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
	public MergedInputNetwork getInputNetwork()
	{
		return input;
	}
	

	/* (non-Javadoc)
	 * @see be.svlandeg.diffany.core.networks.Network#getStringRepresentation()
	 */
	@Override
	public String getStringRepresentation()
	{
		String result = name + ": merged overlapping network calculated from ";
		result += input.getName();
		return result;
	}

}
