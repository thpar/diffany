package be.svlandeg.diffany.examples;

import java.util.Collection;

import org.cytoscape.work.TaskMonitor;

import be.svlandeg.diffany.core.io.EdgeIO;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunConfiguration;
import be.svlandeg.diffany.core.project.RunDiffConfiguration;
import be.svlandeg.diffany.core.project.RunOutput;

/**
 * Generic class for printing an example to the console.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class GenericExample
{
	
	private String name = "";
	protected TaskMonitor taskMonitor;
	
	/**
	 * Allows subclasses to specify a default project for this example. 
	 * 
	 * @return the default project. Returns null if no default was defined.
	 */
	public abstract Project getDefaultProject();
	
	/**
	 * Allows subclasses to specify a default runconfig for this example. 
	 * @param p the project used in this example
	 * 
	 * @return the default configuration. Returns -1 if no default was defined.
	 */
	public abstract int getDefaultRunConfigurationID(Project p);

	/**
	 * Print a network by printing its string representation and its edges.
	 * 
	 * @param n the network to be printed
	 */
	protected void printNetwork(Network n)
	{
		System.out.println(n.getStringRepresentation());
		System.out.println(EdgeIO.writeEdgesToTab(n.getEdges()));
		System.out.println("");
	}

	/**
	 * Print a differential run configuration 
	 */
	protected void printAllNetworks(RunDiffConfiguration rc)
	{
		System.out.println("Reference network : ");
		ReferenceNetwork r = rc.getReferenceNetwork();
		printNetwork(r);

		System.out.println("Condition-specific network(s) : ");
		Collection<ConditionNetwork> cnetworks = rc.getConditionNetworks();
		for (ConditionNetwork c : cnetworks)
		{
			printNetwork(c);
		}
	}

	/**
	 * Print a run configuration
	 */
	protected void printAllNetworks(RunConfiguration rc)
	{
		System.out.println("Input networks : ");

		Collection<InputNetwork> cnetworks = rc.getInputNetworks();
		for (InputNetwork c : cnetworks)
		{
			printNetwork(c);
		}
	}

	/**
	 * Print an entire project.
	 * 
	 * @param p the project to be printed
	 * @param runID the ID of the Run
	 * @param pair whether output pairs should be printed
	 * @param consensusOnly whether only consensus networks should be printed
	 * @param diffOnly whether only differential networks should be printed
	 */
	protected void printAllNetworks(Project p, int runID, boolean pair, boolean consensusOnly, boolean diffOnly)
	{
		RunConfiguration rc = p.getRunConfiguration(runID);
		if (rc.getClass().equals("RunDiffConfiguration"))
		{
			printAllNetworks((RunDiffConfiguration) rc);
		}
		else
		{
			printAllNetworks(rc);
		}
		
		RunOutput output = p.getOutput(runID);
		if (pair)
		{
			System.out.println("Differential network(s) : ");
			for (OutputNetworkPair op : output.getOutputAsPairs())
			{
				printNetwork(op.getDifferentialNetwork());
				printNetwork(op.getConsensusNetwork());
			}
		}
		if (consensusOnly)
		{
			System.out.println("Consensus network(s) only : ");
			for (ConsensusNetwork on : output.getConsensusNetworks())
			{
				printNetwork(on);
			}
		}
		if (diffOnly)
		{
			System.out.println("Differential network(s) only : ");
			for (DifferentialNetwork dn : output.getDifferentialNetworks())
			{
				printNetwork(dn);
			}
		}
	}

	public String getName() 
	{
		return name;
	}

	public void setName(String name) 
	{
		this.name = name;
	}

	
	public void setTaskMonitor(TaskMonitor taskMonitor) 
	{
		this.taskMonitor = taskMonitor;
		
	}
	
	
}
