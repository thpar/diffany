package be.svlandeg.diffany.core.networks;

import java.util.Set;


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
	 * 
	 */
	public InputNetwork(String name, int ID, Set<String> nodeAttributes)
	{
		super(name, ID, nodeAttributes);
	}

	/**
	 * Create a new generic input network.
	 * 
	 * @param networkName the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network
	 * @param nodes the nodes of this network
	 * @param edges the edges of this network
	 */
	public InputNetwork(String networkName, int ID, Set<String> nodeAttributes, Set<Node> nodes, Set<Edge> edges)
	{
		super(networkName, ID, nodeAttributes, nodes, edges);
	}

	@Override
	public String getStringRepresentation()
	{
		return name + ": input network";
	}
}
