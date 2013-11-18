package be.svlandeg.diffany.examples;

import java.util.*;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.*;
import be.svlandeg.diffany.semantics.*;


/** 
 * This class provides examples taken from the Bandyopadhyay et al, Science 2010 paper.

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
		Map<String, Node> nodes = getNodesByNameFigure1C();
		ReferenceNetwork r = getReferenceFigure1C(nodes);
		Set<ConditionNetwork> c = getConditionFigure1C(nodes);
		EdgeOntology eo = new DefaultSymmEdgeOntology();
		NodeMapper nm = new DefaultNodeMapper();
		Project p = new Project(name, r, c, eo, nm);
		return p;
	}

	/**
	 * Get all the nodes neaded for Figure 1C.
	 * @return all nodes in figure 1C, mapped by name.
	 */
	private Map<String, Node> getNodesByNameFigure1C()
	{
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("C", new Node("C"));
		nodes.put("D", new Node("D"));
		nodes.put("E", new Node("E"));
		return nodes;
	}

	/**
	 * Get the reference network depicted in figure 1C.
	 * @return the reference network
	 */
	private ReferenceNetwork getReferenceFigure1C(Map<String, Node> nodes)
	{
		ReferenceNetwork network = new ReferenceNetwork("Untreated");
		network.addEdge(new Edge("negative", nodes.get("A"), nodes.get("D"), true));
		network.addEdge(new Edge("negative", nodes.get("A"), nodes.get("B"), true));
		network.addEdge(new Edge("positive", nodes.get("E"), nodes.get("C"), true));
		return network;
	}

	/**
	 * Get the condition-specific network depicted in figure 1C.
	 * @return the condition-specific network
	 */
	private Set<ConditionNetwork> getConditionFigure1C(Map<String, Node> nodes)
	{
		Set<ConditionNetwork> cnetworks = new HashSet<ConditionNetwork>();

		String description = "treated with MMS";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Treated", conditions);

		network.addEdge(new Edge("negative", nodes.get("A"), nodes.get("D"), true));
		network.addEdge(new Edge("positive", nodes.get("A"), nodes.get("B"), true));
		network.addEdge(new Edge("negative", nodes.get("A"), nodes.get("C"), true));

		cnetworks.add(network);
		return cnetworks;
	}

	/*
	 * Testing the example
	 */
	public static void main(String[] args)
	{
		Bandyopadhyay2010 ex = new Bandyopadhyay2010();
		
		System.out.println("Defining network for Bandyopadhyay2010 figure 1c");
		Project p = ex.getProjectFigure1C();
		
		System.out.println("Calculating differential networks ");
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p);
		
		System.out.println("");
		ex.printAllNetworks(p);
	}
	
}
