package be.svlandeg.diffany.core.networks;

import java.util.Set;

import be.svlandeg.diffany.core.semantics.NodeMapper;


/**
 * A generic input network is not yet defined as a specific Diffany network type, and can at times be used to generically read input data,
 * which still requires pre-processing.
 * 
 * @author Sofie Van Landeghem
 *
 */
public class InputNetwork extends Network
{

	/**
	 * Create a new generic input network.
	 * 
	 * @param name the name of this network
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 * 
	 */
	public InputNetwork(String name, NodeMapper nm)
	{
		super(name, nm);
	}

	/**
	 * Create a new generic input network.
	 * 
	 * @param name the name of this network (should be enforced to be unique within one project)
	 * @param nodes the nodes of this network
	 * @param edges the edges of this network
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 */
	public InputNetwork(String networkName, Set<Node> nodes, Set<Edge> edges, NodeMapper nm)
	{
		super(networkName, nodes, edges, nm);
	}

	@Override
	public String getStringRepresentation()
	{
		return name + ": input network";
	}
}
