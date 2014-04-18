package be.svlandeg.diffany.core.networks.merged;

import java.util.Set;

import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.semantics.NodeMapper;

/**
 * A merged differential network contains differential edges between 2 (or more) networks,
 * one of which is always a 'static' reference network. 
 * It merges the results of one or more normal/small differential networks.
 * 
 * @author Sofie Van Landeghem
 */
public class MergedDifferentialNetwork extends MergedNetwork
{
	
	protected MergedInputNetwork input;
	
	
	/**
	 * Create a new differential network, referring to a merged input network.
	 * 
	 * @param name the name of this network
	 * @param nodes the nodes of this network
	 * @param conditionEdges the edges of this network
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 * @param input the merged input network
	 * 
	 * @throws IllegalArgumentException when the list of condition-specific networks is null or empty,
	 * or when the reference network is null
	 */
	public MergedDifferentialNetwork(String name, Set<Node> nodes, Set<ConditionEdge> conditionEdges, NodeMapper nm, MergedInputNetwork input) 
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
		String result = name + ": merged differential network calculated from ";
		result += input.getName();
		return result;
	}

}
