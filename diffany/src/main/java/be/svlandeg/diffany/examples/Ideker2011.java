package be.svlandeg.diffany.examples;

import java.util.*;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.*;
import be.svlandeg.diffany.semantics.*;


/** 
 * This class provides examples taken from the Ideker et al, Molecular Systems Biology 2011 paper.
 * http://www.nature.com/msb/journal/v8/n1/pdf/msb201199.pdf
 * 
 * @author Sofie Van Landeghem
 */
public class Ideker2011 extends GenericExample
{

	/**
	 * Get a project with the networks as described in figure 3A of this paper.
	 * @return an example project illustrating figure 1C.
	 */
	public Project getProjectFigure3A()
	{
		String name = "Ideker2011_fig3A";
		ReferenceNetwork r = getReferenceFigure3A();
		Set<ConditionNetwork> c = getConditionFigure3A();
		EdgeOntology eo = new DefaultEdgeOntology();
		NodeMapper nm = new DefaultNodeMapper();
		Project p = new Project(name, r, c, eo, nm);
		return p;
	}

	/**
	 * Get the reference network depicted in figure 3A.
	 * @return the reference network
	 */
	private ReferenceNetwork getReferenceFigure3A()
	{
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("D", new Node("D"));
		nodes.put("E", new Node("E"));
		nodes.put("F", new Node("F"));
		
		ReferenceNetwork network = new ReferenceNetwork("Condition 1");
		network.addEdge(new Edge("negative", nodes.get("A"), nodes.get("F"), true));
		network.addEdge(new Edge("negative", nodes.get("A"), nodes.get("D"), true));
		network.addEdge(new Edge("negative", nodes.get("A"), nodes.get("B"), true));
		network.addEdge(new Edge("positive", nodes.get("A"), nodes.get("E"), true));
		return network;
	}

	/**
	 * Get the condition-specific network depicted in figure 3A.
	 * @return the condition-specific network
	 */
	private Set<ConditionNetwork> getConditionFigure3A()
	{
		Set<ConditionNetwork> cnetworks = new HashSet<ConditionNetwork>();

		String description = "Unknown condition";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Condition 2", conditions);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("D", new Node("D"));
		nodes.put("C", new Node("C"));
		nodes.put("F", new Node("F"));

		network.addEdge(new Edge("negative", nodes.get("F"), nodes.get("A"), true));
		network.addEdge(new Edge("positive", nodes.get("B"), nodes.get("A"), true));
		network.addEdge(new Edge("negative", nodes.get("C"), nodes.get("A"), true));
		network.addEdge(new Edge("negative", nodes.get("D"), nodes.get("A"), true));
		
		cnetworks.add(network);
		return cnetworks;
	}

	/**
	 * Testing the example
	 */
	public static void main(String[] args)
	{
		Ideker2011 ex = new Ideker2011();
		
		System.out.println("Defining network for Ideker2011 figure 3A");
		Project p = ex.getProjectFigure3A();
		
		System.out.println("Calculating differential networks ");
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p);
		
		System.out.println("");
		ex.printAllNetworks(p);
	}
}
