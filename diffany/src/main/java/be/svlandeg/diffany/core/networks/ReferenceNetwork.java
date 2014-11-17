package be.svlandeg.diffany.core.networks;

import java.util.Set;


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
	 * 
	 */
	public ReferenceNetwork(String name, int ID, Set<String> nodeAttributes)
	{
		super(name, ID, nodeAttributes);
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
