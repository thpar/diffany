package be.svlandeg.diffany.examples;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.Condition;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.LogEntry;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

/**
 * Testing class that tries to simulate a range of possibilities in activity flow networks
 * and is able to produce differential networks for them.
 * 
 * @author Sofie Van Landeghem
 */
public class ActivityFlowTest extends GenericExample
{
	
	private NodeMapper nm;
	
	/**
	 * Constructor: generates a default {@link NodeMapper} object
	 */
	public ActivityFlowTest()
	{
		nm = new DefaultNodeMapper();
	}
	
	/**
	 * Get a custom project
	 * @return an example project.
	 */
	public Project getTestProject()
	{
		String name = "AF_test";
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo, nm);
		return p;
	}
	
	/**
	 * Add some custom-defined networks to the project.
	 * @param p the input project
	 * @return the resulting configuration ID
	 */
	public int getTestConfiguration(Project p)
	{
		ReferenceNetwork r = getTestReference();
		Set<ConditionNetwork> c = getTestCondition();
		boolean cleanInput = true;
		int ID = p.addRunConfiguration(r, c, cleanInput, null);
		return ID;
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
		nodes.put("G", new Node("G"));
		nodes.put("H", new Node("H"));
		nodes.put("J", new Node("J"));
		nodes.put("K", new Node("K"));
		
		ReferenceNetwork network = new ReferenceNetwork("Condition 1", 1, nm);
		
		// non-alphanumerical chars (punctuation, spaces, ...) should be ignored!
		network.addEdge(new Edge(" positi ve reg ulation", nodes.get("A"), nodes.get("B"), false, 2, false));
		network.addEdge(new Edge("positive regulation", nodes.get("B"), nodes.get("A"), false, 1, false));
		
		network.addEdge(new Edge("negative-regulation", nodes.get("M"), nodes.get("N"), false, 5, false));
		
		network.addEdge(new Edge("positive___regulation", nodes.get("X"), nodes.get("Y"), true, 4, true));
		
		network.addEdge(new Edge("negativeregulation", nodes.get("G"), nodes.get("H"), true, 3, true));
		
		network.addEdge(new Edge("negative_-_regulation", nodes.get("J"), nodes.get("K"), true, 2, false));
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

		ConditionNetwork network = new ConditionNetwork("Condition 2", 2, conditions, nm);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("M", new Node("M"));
		nodes.put("N", new Node("N"));
		nodes.put("S", new Node("S"));
		nodes.put("T", new Node("T"));
		nodes.put("X", new Node("X"));
		nodes.put("Y", new Node("Y"));
		nodes.put("G", new Node("G"));
		nodes.put("H", new Node("H"));
		nodes.put("J", new Node("J"));
		nodes.put("K", new Node("K"));
		
		network.addEdge(new Edge("positive regulation", nodes.get("A"), nodes.get("B"), false, 3, false));
		network.addEdge(new Edge("positive regulation", nodes.get("B"), nodes.get("A"), false, 2, true));
		
		network.addEdge(new Edge("positive regulation", nodes.get("M"), nodes.get("N"), true, 7, false));
		
		network.addEdge(new Edge("positive regulation", nodes.get("S"), nodes.get("T"), false, 1, false));
		
		network.addEdge(new Edge("positive regulation", nodes.get("X"), nodes.get("Y"), false, 2, true));
		
		network.addEdge(new Edge("regulation", nodes.get("H"), nodes.get("G"), true, 7, true));
		
		network.addEdge(new Edge("positive regulation", nodes.get("K"), nodes.get("J"), false, 2, true));
		
		cnetworks.add(network);
		return cnetworks;
	}
	
	/**
	 * Testing the example
	 * @param args the (ignored) input argument list
	 */
	public static void main(String[] args)
	{
		ActivityFlowTest ex = new ActivityFlowTest();
		double cutoff = 0.0;
		
		System.out.println("Defining network for negation test");
		Project p = ex.getTestProject();
		int ID = ex.getTestConfiguration(p);
		
		System.out.println("Calculating differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, cutoff, true, true, 3, true, null);
		
		System.out.println("");
		ex.printAllNetworks(p, ID, true, false, false);
		
		System.out.println("Logs:");
		for (LogEntry l : p.getLogger(ID).getAllLogMessages())
		{
			System.out.println(l);
		}
	}
}
