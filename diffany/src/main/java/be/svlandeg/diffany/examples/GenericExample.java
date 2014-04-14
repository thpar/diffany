package be.svlandeg.diffany.examples;

import java.util.Collection;

import be.svlandeg.diffany.core.io.EdgeIO;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.DifferentialOutput;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunDiffConfiguration;

/**
 * Generic class for printing an example to the console.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class GenericExample
{

	/**
	 * Print a network by printing its string representation and its edges.
	 * @param n the network to be printed
	 */
	private void printNetwork(Network n)
	{
		System.out.println(n.getStringRepresentation());
		System.out.println(EdgeIO.writeEdgesToTab(n.getEdges()));
		System.out.println("");
	}

	/**
	 * Print an entire project.
	 * @param p the project to be printed
	 * @param configurationID the ID of the RunConfiguration
	 */
	protected void printAllNetworks(Project p, int configurationID)
	{
		RunDiffConfiguration rc = (RunDiffConfiguration) p.getRunConfiguration(configurationID);
		
		System.out.println("Reference network : ");
		ReferenceNetwork r = rc.getReferenceNetwork();
		printNetwork(r);

		System.out.println("Condition-specific network(s) : ");
		Collection<ConditionNetwork> cnetworks = rc.getConditionNetworks();
		for (ConditionNetwork c : cnetworks)
		{
			printNetwork(c);
		}
		System.out.println("Differential network(s) : ");
		Collection<DifferentialOutput> dnetworks = rc.getDifferentialOutputs();
		
		for (DifferentialOutput d : dnetworks)
		{
			OutputNetworkPair pair = d.getOutputAsPair();
			
			printNetwork(pair.getDifferentialNetwork());
			printNetwork(pair.getOverlappingNetwork());
		}
	}
	
	/**
	 * Print all overlap networks in an entire project.
	 * @param p the project to be printed
	 * @param configurationID the ID of the RunConfiguration
	 */
	protected void printAllOverlapNetworks(Project p, int configurationID)
	{
		RunDiffConfiguration rc = (RunDiffConfiguration) p.getRunConfiguration(configurationID);
		
		System.out.println("Reference network : ");
		ReferenceNetwork r = rc.getReferenceNetwork();
		printNetwork(r);

		System.out.println("Condition-specific network(s) : ");
		Collection<ConditionNetwork> cnetworks = rc.getConditionNetworks();
		for (ConditionNetwork c : cnetworks)
		{
			printNetwork(c);
		}
		System.out.println("Overlap network(s) : ");
		Collection<DifferentialOutput> dnetworks = rc.getDifferentialOutputs();
		
		for (DifferentialOutput d : dnetworks)
		{
			printNetwork(d.getOverlappingNetwork());
		}
	}

}
