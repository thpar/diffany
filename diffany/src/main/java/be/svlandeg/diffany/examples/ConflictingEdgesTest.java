package be.svlandeg.diffany.examples;

import java.util.*;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.*;
import be.svlandeg.diffany.semantics.*;

/**
 * Testing class that simulates input networks with edge conflicts, 
 * e.g. positive and negative regulation at the same time.
 * 
 * @author Sofie Van Landeghem
 */
public class ConflictingEdgesTest extends GenericExample
{
	
	/**
	 * Get a project with some custom-defined networks.
	 * @return an example project.
	 */
	public Project getTestProject()
	{
		String name = "Conflict_test";
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
		nodes.put("G", new Node("G"));
		nodes.put("H", new Node("H"));
		nodes.put("J", new Node("K"));
		nodes.put("K", new Node("J"));
		
		ReferenceNetwork network = new ReferenceNetwork("Reference");
		
		network.addEdge(new Edge("positive regulation", nodes.get("A"), nodes.get("B"), false, 2, false));
		network.addEdge(new Edge("ptm", nodes.get("A"), nodes.get("B"), false, 5, false));
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

		ConditionNetwork network = new ConditionNetwork("Condition-specific network", conditions);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("G", new Node("G"));
		nodes.put("H", new Node("H"));
		nodes.put("J", new Node("K"));
		nodes.put("K", new Node("J"));
		
		network.addEdge(new Edge("positive regulation", nodes.get("A"), nodes.get("B"), false, 8, false));
		
		cnetworks.add(network);
		return cnetworks;
	}
	
	/**
	 * Testing the example
	 */
	public static void main(String[] args)
	{
		ConflictingEdgesTest ex = new ConflictingEdgesTest();
		double cutoff = 0.0;
		
		System.out.println("Defining network for conflict test");
		Project p = ex.getTestProject();
		
		System.out.println("Calculating differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, cutoff);
		
		System.out.println("");
		ex.printAllNetworks(p);
	}

}
