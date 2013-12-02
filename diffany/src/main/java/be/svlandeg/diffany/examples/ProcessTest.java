package be.svlandeg.diffany.examples;

import java.util.*;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.*;
import be.svlandeg.diffany.semantics.*;

/**
 * Testing class that tries to simulate a range of possibilities in process networks
 * and is able to produce differential networks for them.
 * 
 * @author Sofie Van Landeghem
 */
public class ProcessTest extends GenericExample
{
	
	/**
	 * Get a project with some custom-defined networks.
	 * @return an example project.
	 */
	public Project getTestProject()
	{
		String name = "Process_test";
		ReferenceNetwork r = getTestReference();
		Set<ConditionNetwork> c = getTestCondition();
		EdgeOntology eo = new DefaultEdgeOntology();
		NodeMapper nm = new DefaultNodeMapper();
		Project p = new Project(name, r, c, eo, nm);
		return p;
	}
	
	/**
	 * Get a custom-defined reference network.
	 * @return the reference network
	 */
	private ReferenceNetwork getTestReference()
	{
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("M", new Node("M"));
		nodes.put("N", new Node("N"));
		nodes.put("S", new Node("S"));
		nodes.put("T", new Node("T"));
		nodes.put("X", new Node("X"));
		nodes.put("Y", new Node("Y"));
		
		ReferenceNetwork network = new ReferenceNetwork("Condition 1");
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 2, false));
		network.addEdge(new Edge("ppi", nodes.get("S"), nodes.get("T"), true, 3, false));
		network.addEdge(new Edge("ppi", nodes.get("X"), nodes.get("Y"), true, 4, false));
		return network;
	}
	
	/**
	 * Get the custom-defined condition-specific network.
	 * @return the condition-specific network
	 */
	private Set<ConditionNetwork> getTestCondition()
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
		nodes.put("M", new Node("M"));
		nodes.put("N", new Node("N"));
		nodes.put("S", new Node("S"));
		nodes.put("T", new Node("T"));
		nodes.put("X", new Node("X"));
		nodes.put("Y", new Node("Y"));

		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 3, false));
		network.addEdge(new Edge("ppi", nodes.get("M"), nodes.get("N"), true, 3, false));
		network.addEdge(new Edge("ppi", nodes.get("S"), nodes.get("T"), true, 2, false));
		
		cnetworks.add(network);
		return cnetworks;
	}
	
	/**
	 * Testing the example
	 */
	public static void main(String[] args)
	{
		ProcessTest ex = new ProcessTest();
		double cutoff = 0.0;
		
		System.out.println("Defining network for weight test");
		Project p = ex.getTestProject();
		
		System.out.println("Calculating differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, cutoff);
		
		System.out.println("");
		ex.printAllNetworks(p);
	}

}
