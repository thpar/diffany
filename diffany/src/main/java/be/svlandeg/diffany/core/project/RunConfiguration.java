package be.svlandeg.diffany.core.project;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.core.networks.InputNetwork;


/**
 * A RunConfiguration defines the necessary networks needed as input for the Diffany algorithms,
 * and should always be used in the context of a bigger {@link Project}.
 * 
 * This configuration can be used to calculate overlap networks, but differential networks can only be calculated with a RunDiffConfiguration object.
 * 
 * @author Sofie Van Landeghem
 */
public class RunConfiguration
{
	
	protected Set<InputNetwork> inputNetworks;
	protected Set<DifferentialOutput> differentialOutputs;
	
	
	/**
	 * Create a new configuration with a set of input networks. 
	 * The output result set is initialized to be empty.
	 * 
	 * @param inputNetworks the input networks (not null or empty!)
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public RunConfiguration(Set<InputNetwork> inputNetworks)
	{
		setInputs(inputNetworks);
		cleanOutputResults();
	}
	
	/**
	 * Set the input networks in this configuration.
	 * (currently not a public method - changes to it would influence the differential networks (TODO v3.0))
	 * 
	 * @param inputNetworks the input networks (not null or empty!)
	 * @throws IllegalArgumentException if the set is null or empty
	 */
	protected void setInputs(Set<InputNetwork> inputNetworks) throws IllegalArgumentException
	{
		if (inputNetworks == null || inputNetworks.isEmpty())
		{
			String errormsg = "The set of input networks should not be null or empty!";
			throw new IllegalArgumentException(errormsg);
		}
		this.inputNetworks = inputNetworks;
	}
	
	/**
	 * Add a differential output network to this configuration, with or without prior cleaning of the previous result set.
	 * 
	 * @param differential a new differential network
	 * @param clean whether or not to first clean the existing result set
	 */
	public void addOutputResult(DifferentialOutput output, boolean clean)
	{
		if (clean)
		{
			cleanOutputResults();
		}
		differentialOutputs.add(output);
	}
	
	/**
	 * Remove all previous differential output networks from this configuration.
	 */
	public void cleanOutputResults()
	{
		differentialOutputs = new HashSet<DifferentialOutput>();
	}
	
	
	/**
	 * Get all input network(s): 1 or many. 
	 * 
	 * @return the input networks in this configuration (1 or more, never null or empty))
	 */
	public Collection<InputNetwork> getInputNetworks()
	{
		return inputNetworks;
	}

	/**
	 * Get the output networks in the configuration: 0, 1 or more
	 * 
	 * @return the differential output networks in this configuration (if any, otherwise empty set, but never null)
	 */
	public Collection<DifferentialOutput> getDifferentialOutputs()
	{
		return differentialOutputs;
	}

}
