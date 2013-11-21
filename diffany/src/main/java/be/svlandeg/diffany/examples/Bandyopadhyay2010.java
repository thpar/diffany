package be.svlandeg.diffany.examples;

import java.util.*;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.*;
import be.svlandeg.diffany.semantics.*;

/** 
 * This class provides examples taken from the Bandyopadhyay et al, Science 2010 paper.
 * http://www.sciencemag.org/content/330/6009/1385.full.pdf
 * 
 * @author Sofie Van Landeghem
 */
public class Bandyopadhyay2010 extends GenericExample
{

	/**
	 * Get a project with the networks as described in figure 1C of this paper.
	 * @return an example project illustrating figure 1C.
	 */
	public Project getProjectFigure1C()
	{
		String name = "Bandyopadhyay2010_fig1C";
		ReferenceNetwork r = getReferenceFigure1C();
		Set<ConditionNetwork> c = getConditionFigure1C();
		EdgeOntology eo = new DefaultEdgeOntology();
		NodeMapper nm = new DefaultNodeMapper();
		Project p = new Project(name, r, c, eo, nm);
		return p;
	}


	/**
	 * Get the reference network depicted in figure 1C.
	 * @return the reference network
	 */
	private ReferenceNetwork getReferenceFigure1C()
	{
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("C", new Node("C"));
		nodes.put("D", new Node("D"));
		nodes.put("E", new Node("E"));
		
		ReferenceNetwork network = new ReferenceNetwork("Untreated");
		network.addEdge(new Edge("negative", nodes.get("A"), nodes.get("D"), true, 1.1, false));
		network.addEdge(new Edge("negative", nodes.get("A"), nodes.get("B"), true, 0.3, false));
		network.addEdge(new Edge("positive", nodes.get("E"), nodes.get("C"), true, 0.8, false));
		return network;
	}

	/**
	 * Get the condition-specific network depicted in figure 1C.
	 * @return the condition-specific network
	 */
	private Set<ConditionNetwork> getConditionFigure1C()
	{
		Set<ConditionNetwork> cnetworks = new HashSet<ConditionNetwork>();

		String description = "treated with MMS";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Treated", conditions);

		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("C", new Node("C"));
		nodes.put("D", new Node("D"));
		
		network.addEdge(new Edge("negative", nodes.get("A"), nodes.get("D"), true, 0.9, false));
		network.addEdge(new Edge("positive", nodes.get("A"), nodes.get("B"), true, 0.4, false));
		network.addEdge(new Edge("negative", nodes.get("A"), nodes.get("C"), true, 0.7, false));

		cnetworks.add(network);
		return cnetworks;
	}

	/**
	 * Testing the example
	 */
	public static void main(String[] args)
	{
		Bandyopadhyay2010 ex = new Bandyopadhyay2010();
		double cutoff = 0.25;
		
		System.out.println("Defining network for Bandyopadhyay2010 figure 1c");
		Project p = ex.getProjectFigure1C();
		
		System.out.println("Calculating differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, cutoff);
		
		System.out.println("");
		ex.printAllNetworks(p);
	}
	
}
