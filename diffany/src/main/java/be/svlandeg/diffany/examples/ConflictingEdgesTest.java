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

	private NodeMapper nm;
	
	/**
	 * Constructor: generates a default {@link NodeMapper} object
	 */
	public ConflictingEdgesTest()
	{
		nm = new DefaultNodeMapper();
	}
	
	/**
	 * Get a custom project.
	 * @return an example project.
	 */
	public Project getTestProject()
	{
		String name = "Conflict_test";
		EdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo, nm);
		return p;
	}
	
	/**
	 * Add some custom-defined networks to the project.
	 * @return the resulting configuration ID.
	 */
	public int getTestConfiguration(Project p)
	{
		ReferenceNetwork r = getTestReference();
		Set<ConditionNetwork> c = getTestCondition();
		RunConfiguration rc = new RunConfiguration (r, c);
		int ID = p.addRunConfiguration(rc);
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
		nodes.put("G", new Node("G"));
		nodes.put("H", new Node("H"));
		nodes.put("J", new Node("J"));
		nodes.put("K", new Node("K"));
		
		ReferenceNetwork network = new ReferenceNetwork("Reference", nm);
		
		network.addEdge(new Edge("positive regulation", nodes.get("A"), nodes.get("B"), false, 2, false));
		network.addEdge(new Edge("ptm", nodes.get("A"), nodes.get("B"), false, 5, false));
		network.addEdge(new Edge("somerandomInteraction", nodes.get("A"), nodes.get("B"), false, 4, false));
		
		network.addEdge(new Edge("ptm", nodes.get("G"), nodes.get("H"), false, 3, false));
		network.addEdge(new Edge("phosphorylation", nodes.get("G"), nodes.get("H"), false, 1, false));
		network.addEdge(new Edge("regulation", nodes.get("G"), nodes.get("H"), false, 7, false));
		
		network.addEdge(new Edge("ptm", nodes.get("J"), nodes.get("K"), false, 4, true));
		network.addEdge(new Edge("ubiquitination", nodes.get("J"), nodes.get("K"), true, 2, false));
		network.addEdge(new Edge("regulation", nodes.get("J"), nodes.get("K"), false, 7, false));
		network.addEdge(new Edge("regulation", nodes.get("J"), nodes.get("K"), false, 5, true));
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

		ConditionNetwork network = new ConditionNetwork("Condition-specific network", conditions, nm);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("G", new Node("G"));
		nodes.put("H", new Node("H"));
		nodes.put("J", new Node("J"));
		nodes.put("K", new Node("K"));
		
		network.addEdge(new Edge("positive regulation", nodes.get("A"), nodes.get("B"), false, 8, false));
		
		network.addEdge(new Edge("phosphorylation", nodes.get("G"), nodes.get("H"), false, 2, false));
		network.addEdge(new Edge("catalysis", nodes.get("G"), nodes.get("H"), false, 4, false));
		
		network.addEdge(new Edge("ptm", nodes.get("J"), nodes.get("K"), true, 3, false));
		network.addEdge(new Edge("ptm", nodes.get("J"), nodes.get("K"), false, 6, false));
		
		network.addEdge(new Edge("positive_regulation", nodes.get("J"), nodes.get("K"), false, 3, true));
		network.addEdge(new Edge("inhibition", nodes.get("J"), nodes.get("K"), true, 4, false));
		
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
		int ID = ex.getTestConfiguration(p);
		
		System.out.println("Calculating differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, cutoff);
		
		System.out.println("");
		ex.printAllNetworks(p, ID);
		
		System.out.println("Logs:");
		for (String l : p.getLogger(ID).getAllLogMessages())
		{
			System.out.println(l);
		}
	}

}
