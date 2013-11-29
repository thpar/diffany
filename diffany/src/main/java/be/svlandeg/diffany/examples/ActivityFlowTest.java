package be.svlandeg.diffany.examples;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.Condition;
import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.Node;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * Testing class that tries to simulate a range of possibilities in activity flow networks
 * and is able to produce differential networks for them.
 * 
 * @author Sofie Van Landeghem
 */
public class ActivityFlowTest extends GenericExample
{
	/**
	 * Get a project with some custom-defined networks.
	 * @return an example project.
	 */
	public Project getTestProject()
	{
		String name = "Negation_test";
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
		network.addEdge(new Edge("positive regulation", nodes.get("A"), nodes.get("B"), false, 2, false));
		network.addEdge(new Edge("positive regulation", nodes.get("B"), nodes.get("A"), false, 1, false));
		
		network.addEdge(new Edge("negative regulation", nodes.get("M"), nodes.get("N"), false, 5, false));
		
		network.addEdge(new Edge("positive regulation", nodes.get("S"), nodes.get("T"), false, 5, true));
		
		network.addEdge(new Edge("positive regulation", nodes.get("X"), nodes.get("Y"), true, 4, true));
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

		network.addEdge(new Edge("positive regulation", nodes.get("A"), nodes.get("B"), false, 3, false));
		network.addEdge(new Edge("positive regulation", nodes.get("B"), nodes.get("A"), false, 2, true));
		
		network.addEdge(new Edge("positive regulation", nodes.get("M"), nodes.get("N"), true, 7, false));
		
		network.addEdge(new Edge("positive regulation", nodes.get("S"), nodes.get("T"), false, 1, true));
		
		network.addEdge(new Edge("positive regulation", nodes.get("X"), nodes.get("Y"), false, 2, true));
		
		cnetworks.add(network);
		return cnetworks;
	}
	
	/**
	 * Testing the example
	 */
	public static void main(String[] args)
	{
		ActivityFlowTest ex = new ActivityFlowTest();
		double cutoff = 0.0;
		
		System.out.println("Defining network for negation test");
		Project p = ex.getTestProject();
		
		System.out.println("Calculating differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, cutoff);
		
		System.out.println("");
		ex.printAllNetworks(p);
	}
}
