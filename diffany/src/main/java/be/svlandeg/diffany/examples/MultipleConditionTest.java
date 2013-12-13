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
 * Testing class that tries to simulate a use-case of calculating a differential network
 * from one reference network and multiple condition-specific networks.
 * 
 * @author Sofie Van Landeghem
 */
public class MultipleConditionTest extends GenericExample
{
	/**
	 * Get a project with some custom-defined networks.
	 * @return an example project.
	 */
	public Project getTestProject()
	{
		String name = "Multiple_test";
		ReferenceNetwork r = getTestReference();
		Set<ConditionNetwork> c = getTestConditions();
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
		nodes.put("C", new Node("C"));
		nodes.put("W", new Node("W"));
		nodes.put("Z", new Node("Z"));
		
		ReferenceNetwork network = new ReferenceNetwork("Reference");
		
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 2, false));
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("C"), true, 2, false));
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("Z"), true, 2, false));
		network.addEdge(new Edge("ppi", nodes.get("Z"), nodes.get("W"), true, 2, false));
		
		nodes.put("M", new Node("M"));
		nodes.put("N", new Node("N"));
		network.addEdge(new Edge("phosphorylation", nodes.get("M"), nodes.get("N"), true, 2, false));

		return network;
	}
	
	/**
	 * Get the custom-defined condition-specific network.
	 * @return the condition-specific network
	 */
	private Set<ConditionNetwork> getTestConditions()
	{
		Set<ConditionNetwork> cnetworks = new HashSet<ConditionNetwork>();
		cnetworks.add(getFirstCondition());
		cnetworks.add(getSecondCondition());
		return cnetworks;
	}
	
	/**
	 * Get the custom-defined condition-specific network.
	 * @return the condition-specific network
	 */
	private ConditionNetwork getFirstCondition()
	{
		String description = "Salt stress";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Condition 1", conditions);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("C", new Node("C"));
		nodes.put("D", new Node("D"));
		nodes.put("F", new Node("F"));
		nodes.put("W", new Node("W"));
		nodes.put("Z", new Node("Z"));

		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 2, false));
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("C"), true, 2, false));
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("Z"), true, 2, false));
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("D"), true, 2, false));
		network.addEdge(new Edge("ppi", nodes.get("D"), nodes.get("F"), true, 2, false));
		
		nodes.put("M", new Node("M"));
		nodes.put("N", new Node("N"));
		network.addEdge(new Edge("phosphorylation", nodes.get("M"), nodes.get("N"), true, 6, false));

		return network;
	}
	
	/**
	 * Get the custom-defined condition-specific network.
	 * @return the condition-specific network
	 */
	private ConditionNetwork getSecondCondition()
	{
		String description = "Draught stress";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Condition 2", conditions);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("C", new Node("C"));
		nodes.put("D", new Node("D"));
		nodes.put("E", new Node("E"));
		nodes.put("W", new Node("W"));
		nodes.put("Z", new Node("Z"));

		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 2, false));
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("C"), true, 2, false));
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("D"), true, 2, false));
		network.addEdge(new Edge("ppi", nodes.get("D"), nodes.get("E"), true, 2, false));
				
		nodes.put("M", new Node("M"));
		nodes.put("N", new Node("N"));
		network.addEdge(new Edge("phosphorylation", nodes.get("M"), nodes.get("N"), true, 8, false));

		return network;
	}
	
	/**
	 * Testing the example
	 */
	public static void main(String[] args)
	{
		MultipleConditionTest ex = new MultipleConditionTest();
		double cutoff = 0.0;
		
		System.out.println("Defining network for weight test");
		Project p = ex.getTestProject();
		
		System.out.println("Calculating differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateOneDifferentialNetwork(p, cutoff);
		
		System.out.println("");
		ex.printAllNetworks(p);
	}
}
