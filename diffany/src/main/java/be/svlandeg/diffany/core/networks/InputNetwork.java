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
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network - can be left empty or null
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 * 
	 */
	public InputNetwork(String name, int ID, Set<String> nodeAttributes, NodeMapper nm)
	{
		super(name, ID, nodeAttributes, nm);
	}

	/**
	 * Create a new generic input network.
	 * 
	 * @param networkName the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network
	 * @param nodes the nodes of this network
	 * @param edges the edges of this network
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 */
	public InputNetwork(String networkName, int ID, Set<String> nodeAttributes, Set<Node> nodes, Set<Edge> edges, NodeMapper nm)
	{
		super(networkName, ID, nodeAttributes, nodes, edges, nm);
	}

	@Override
	public String getStringRepresentation()
	{
		return name + ": input network";
	}
}
