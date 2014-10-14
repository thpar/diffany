package be.svlandeg.diffany.core.networks;

import java.util.Set;

import be.svlandeg.diffany.core.semantics.NodeMapper;


/**
 * A reference network is used as comparison against condition-dependent networks.
 * It can be seen as a 'static' network with unspecified/unknown/wild-type conditions.
 * 
 * @author Sofie Van Landeghem
 *
 */
public class ReferenceNetwork extends InputNetwork
{
	
	/**
	 * Create a new static reference network.
	 * 
	 * @param name the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network - can be left empty or null
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 * 
	 */
	public ReferenceNetwork(String name, int ID, Set<String> nodeAttributes, NodeMapper nm)
	{
		super(name, ID, nodeAttributes, nm);
	}

	/* (non-Javadoc)
	 * @see be.svlandeg.diffany.core.networks.Network#getStringRepresentation()
	 */
	@Override
	public String getStringRepresentation()
	{
		return name + ": input reference network";
	}

	

	
}
