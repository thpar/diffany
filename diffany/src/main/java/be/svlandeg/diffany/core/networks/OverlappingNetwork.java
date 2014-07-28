package be.svlandeg.diffany.core.networks;

import java.util.Set;

import be.svlandeg.diffany.core.semantics.NodeMapper;


/**
 * An overlapping network is the counterpart of a differential network: containing everything but the differential edges 
 * between 2 (or more) networks, an overlapping network stores the overlap or common edges between them.
 * Such a network may be useful to study 'house-keeping' interactions.
 * 
 * In contrast to a differential network, an overlappingnetwork can also be between two condition-specific networks 
 * and doesn't necessarily need a reference network.
 * 
 * @author Sofie Van Landeghem
 */
public class OverlappingNetwork extends Network
{
	
	protected Set<Network> originalNetworks;

	/**
	 * Create a new overlapping network, referring to the original networks it was created from.
	 * 
	 * @param name the name of this network
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param originalNetworks the original networks (at least 2!)
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 * 
	 * @throws IllegalArgumentException when the list of original networks is null or contains less than 2 networks
	 */
	public OverlappingNetwork(String name, int ID, Set<Network> originalNetworks, NodeMapper nm) throws IllegalArgumentException
	{
		super(name, ID, nm);
		if (originalNetworks == null || originalNetworks.size() < 2)
		{
			String errormsg = "Please define at least 2 original networks!";
			throw new IllegalArgumentException(errormsg);
		}
		this.originalNetworks = originalNetworks;
	}
	
	@Override
	public String getStringRepresentation()
	{
		String result = name + ": overlapping network between ";
		for (Network n : originalNetworks)
		{
			result += n.getName() + " / ";
		}
		result = result.substring(0, result.length() - 2);
		return result;
	}

}
