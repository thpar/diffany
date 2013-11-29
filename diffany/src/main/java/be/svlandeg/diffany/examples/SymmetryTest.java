package be.svlandeg.diffany.examples;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.*;
import be.svlandeg.diffany.semantics.*;

/**
 * Toy testing class to experiment with edge symmetry in custom networks.
 * 
 * @author Sofie Van Landeghem
 */
public class SymmetryTest extends GenericExample
{
	/**
	 * Get a project with some custom-defined networks.
	 * @return an example project.
	 */
	public Project getTestProject()
	{
		String name = "Symmetry_test";
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
		network.addEdge(new Edge("positive regulation", nodes.get("A"), nodes.get("B"), false, 2));
		network.addEdge(new Edge("positive regulation", nodes.get("B"), nodes.get("A"), false, 1));
		
		network.addEdge(new Edge("positive regulation", nodes.get("M"), nodes.get("N"), true, 5));
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

		network.addEdge(new Edge("positive regulation", nodes.get("B"), nodes.get("A"), false, 3));
		network.addEdge(new Edge("negative regulation", nodes.get("M"), nodes.get("N"), true, 4));
		
		cnetworks.add(network);
		return cnetworks;
	}
	
	/**
	 * Testing the example
	 */
	public static void main(String[] args)
	{
		SymmetryTest ex = new SymmetryTest();
		double cutoff = 0.0;
		
		System.out.println("Defining network for symmetry test");
		Project p = ex.getTestProject();
		
		System.out.println("Calculating differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, cutoff);
		
		System.out.println("");
		ex.printAllNetworks(p);
	}

}
