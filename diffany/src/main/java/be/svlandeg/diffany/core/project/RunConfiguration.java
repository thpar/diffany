package be.svlandeg.diffany.core.project;

import java.util.Collection;
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
	protected DifferentialOutput differentialOutput;
	
	protected int overlapNo_cutoff;
	
	
	/**
	 * Create a new configuration with a set of input networks. The required overlap cutoff is by default set to the size of this set.
	 * The output result set is initialized to be empty.
	 * 
	 * @param inputNetworks the input networks (not null or empty!)
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public RunConfiguration(Set<InputNetwork> inputNetworks)
	{
		setInputs(inputNetworks);
		overlapNo_cutoff = inputNetworks.size();
		differentialOutput = new DifferentialOutput(); 
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
	 * Alter the required number of overlapping edges needed before the edge will be present in the overlap network.
	 * @param overlapNo_cutoff
	 */
	public void setOverlapCutoff(int overlapNo_cutoff)
	{
		this.overlapNo_cutoff = overlapNo_cutoff;
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
	 * Get the output produced for this configuration.
	 * 
	 * @return the differential output networks in this configuration
	 */
	public DifferentialOutput getDifferentialOutput()
	{
		return differentialOutput;
	}
	
	/**
	 * Return the required number of overlapping edges needed before the edge will be present in the overlap network.
	 * @return the minimal number of required overlapping edges, which is the size of the input set unless specifically stated otherwise
	 */
	public int getOverlapCutoff()
	{
		return overlapNo_cutoff;
	}

}
