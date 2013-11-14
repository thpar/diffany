package be.svlandeg.diffany.concepts;

import java.util.Set;


/**
 * A shared network is the counterpart of a differential network: containing everything but the differential edges between 2 (or more) networks,
 * a shared network stores the overlap or common edges between them.
 * Such a network may be useful to study 'house-keeping' interactions.
 * 
 * @author Sofie Van Landeghem
 */
public class SharedNetwork extends Network
{
	
	protected Set<Network> originalNetworks;

	/**
	 * Create a new shared network, referring to the original networks it was created from.
	 * 
	 * @param name the name of this network
	 * @param originalNetworks the original networks (at least 2!)
	 * @throws IllegalArgumentException when the list of original networks is null or contains less than 2 networks
	 */
	public SharedNetwork(String name, Set<Network> originalNetworks) throws IllegalArgumentException
	{
		super(name);
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
		String result = name + ": shared network between ";
		for (Network n : originalNetworks)
		{
			result += n.getName() + " / ";
		}
		result = result.substring(0, result.length() - 2);
		return result;
	}

}
