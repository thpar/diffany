package be.svlandeg.diffany.core.networks;

import be.svlandeg.diffany.core.semantics.NodeMapper;


/**
 * A reference network is used as comparison against condition-dependent networks.
 * It can be seen as a 'static' network with unspecified/unknown/wild-type conditions.
 * 
 * @author Sofie Van Landeghem
 *
 */
public class ReferenceNetwork extends Network
{
	
	/**
	 * Create a new static reference network.
	 * 
	 * @param name the name of this network
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 * 
	 */
	public ReferenceNetwork(String name, NodeMapper nm)
	{
		super(name, nm);
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
