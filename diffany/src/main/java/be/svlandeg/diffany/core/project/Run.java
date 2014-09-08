package be.svlandeg.diffany.core.project;

import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.core.networks.InputNetwork;


/**
 * A run consists of a {@link RunConfiguration} which defines the networks that can be used as input for the Diffany algorithms. 
 * 
 * Further, there is a RunOutput and a Logger object.
 * 
 * @author Sloffie
 */
public class Run
{
	
	protected Project p;
	protected int runID;
	
	protected RunConfiguration configuration;
	protected RunOutput output;
	protected Boolean type;	// if true: can do differential
	protected Logger logger;
	
	/**
	 * Create a new run, belonging to a specific project.
	 * 
	 * @param p the project
	 * @param runID the ID of this run within the project
	 * @param configuration the configuration of the input networks
	 * @param type whether or not this run can be used for the calculation of differential networks
	 * @param logger the logger object which stores important messages
	 */
	public Run(Project p, int runID, RunConfiguration configuration, Boolean type, Logger logger)
	{
		this.p = p;
		this.runID = runID;
		this.configuration = configuration;
		this.type = type;
		this.logger = logger;
		output = new RunOutput(p, runID);
	}
	
	/**
	 * Private method that checks whether there are no conflicting network IDs in this run.
	 * @return whether or not all IDs of the input and output networks in this run are unique
	 */
	protected boolean checkInputIDs()
	{
		boolean allOK = true;
		
		Set<Integer> readIDs = new HashSet<Integer>();
		for (InputNetwork input : configuration.inputNetworks)
		{
			int ID = input.getID();
			if (readIDs.contains(ID))
			{
				allOK = false;
			}
			readIDs.add(ID);
		}
		return allOK;
	}
	
	/**
	 * Private method that checks whether a specific ID of an output network can be used in this run
	 * @param outputID the proposed output network ID
	 * @return whether or not the provided output network ID can be used for this run
	 */
	protected boolean checkoutputID(int outputID)
	{
		Set<Integer> readIDs = new HashSet<Integer>();
		for (InputNetwork input : configuration.inputNetworks)
		{
			int ID = input.getID();
			readIDs.add(ID);
		}
		return (! readIDs.contains(outputID));
	}

}
