package be.svlandeg.diffany.examples;

import java.util.Collection;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.concepts.RunConfiguration;
import be.svlandeg.diffany.io.EdgeIO;

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
		RunConfiguration rc = p.getRunConfiguration(configurationID);
		
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
		Collection<DifferentialNetwork> dnetworks = rc.getDifferentialNetworks();
		for (DifferentialNetwork d : dnetworks)
		{
			printNetwork(d);
			printNetwork(d.getOverlappingNetwork());
		}
	}

}
