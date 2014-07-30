package be.svlandeg.diffany.examples;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.Condition;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.LogEntry;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

/** 
 * This class provides examples to benchmark the fuzzy overlap functionality.  
 * 
 * @author Sofie Van Landeghem
 */
public class FuzzyOverlap extends GenericExample
{

	private NodeMapper nm;
	
	/**
	 * Constructor: generates a default {@link NodeMapper} object
	 */
	public FuzzyOverlap()
	{
		nm = new DefaultNodeMapper();
	}
	
	/**
	 * Get a custom project.
	 * @return an example project illustrating figure 1C.
	 */
	public Project getProject()
	{
		String name = "FuzzyOverlap";
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo, nm);
		return p;
	}
	
	/**
	 * Add some custom-defined networks to the project: 1 reference network and 3 condition-specific.
	 * @param p the fuzzy project
	 * @param overlapCutoff the cutoff for overlap %
	 * @return the resulting configuration ID.
	 */
	public int getTestConfigurationWithReference(Project p, int overlapCutoff)
	{
		ReferenceNetwork r = getReference();
		Set<ConditionNetwork> c = new HashSet<ConditionNetwork>();
		c.add(getCondition1());
		c.add(getCondition2());
		c.add(getCondition3());	
		int ID = p.addRunConfiguration(r, c, overlapCutoff);
		return ID;
	}
	
	/**
	 * Add some custom-defined networks to the project: 4 condition-specific networks, no reference.
	 * @param p the fuzzy project
	 * @param overlapCutoff the cutoff for overlap %
	 * @return the resulting configuration ID.
	 */
	public int getTestConfigurationWithoutReference(Project p, int overlapCutoff)
	{
		Set<InputNetwork> c = new HashSet<InputNetwork>();
		c.add(getCondition0());	
		c.add(getCondition1());	
		c.add(getCondition2());
		c.add(getCondition3());	
		int ID = p.addRunConfiguration(c, overlapCutoff);
		return ID;
	}


	/**
	 * Get the reference network useful for testing
	 * @return the reference network
	 */
	private ReferenceNetwork getReference()
	{
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		
		nodes.put("X", new Node("X"));
		nodes.put("Y", new Node("Y"));
		
		nodes.put("M", new Node("M"));
		nodes.put("N", new Node("N"));
		
		ReferenceNetwork network = new ReferenceNetwork("Fuzzy reference network", 1, nm);
		
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 0.8));
		
		network.addEdge(new Edge("negative regulation", nodes.get("X"), nodes.get("Y"), false, 0.5));
		
		network.addEdge(new Edge("ptm", nodes.get("M"), nodes.get("N"), false, 0.5, true));
		
		return network;
	}
	
	/**
	 * Get the first condition-specific network
	 * @return the first condition-specific network
	 */
	private ConditionNetwork getCondition1()
	{
		String description = "treated with MMS 1";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Condition network 1", 11, conditions, nm);

		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		
		nodes.put("X", new Node("X"));
		nodes.put("Y", new Node("Y"));
		
		nodes.put("M", new Node("M"));
		nodes.put("N", new Node("N"));
		
		network.addEdge(new Edge("colocalization", nodes.get("A"), nodes.get("B"), false, 0.6));
		
		network.addEdge(new Edge("positive regulation", nodes.get("X"), nodes.get("Y"), false, 0.8));
		
		network.addEdge(new Edge("phosphorylation", nodes.get("M"), nodes.get("N"), false, 0.3, true));

		return network;
	}

	/**
	 * Get the second condition-specific network
	 * @return the second condition-specific network
	 */
	private ConditionNetwork getCondition2()
	{
		String description = "treated with MMS 2";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Condition network 2", 12, conditions, nm);

		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		
		nodes.put("X", new Node("X"));
		nodes.put("Y", new Node("Y"));
		
		nodes.put("M", new Node("M"));
		nodes.put("N", new Node("N"));
		
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), false, 0.3));
		
		network.addEdge(new Edge("regulation", nodes.get("X"), nodes.get("Y"), false, 0.6));
		
		network.addEdge(new Edge("ptm", nodes.get("M"), nodes.get("N"), false, 0.6, true));

		return network;
	}
	
	/**
	 * Get the third condition-specific network
	 * @return the third condition-specific network
	 */
	private ConditionNetwork getCondition3()
	{
		String description = "treated with MMS 3";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Condition network 3", 13, conditions, nm);

		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		
		nodes.put("X", new Node("X"));
		nodes.put("Y", new Node("Y"));
		
		nodes.put("M", new Node("M"));
		nodes.put("N", new Node("N"));
		
		network.addEdge(new Edge("ppi", nodes.get("B"), nodes.get("A"), false, 0.4));
		
		network.addEdge(new Edge("positive regulation", nodes.get("X"), nodes.get("Y"), false, 0.3));
		
		network.addEdge(new Edge("phosphorylation", nodes.get("M"), nodes.get("N"), false, 0.7, false));

		return network;
	}
	
	/**
	 * Get the "zeroeth" condition-specific network - to be used when there is no reference network
	 * @return the "zeroeth" condition-specific network
	 */
	private ConditionNetwork getCondition0()
	{
		String description = "treated with MMS 0";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Condition network 0", 10, conditions, nm);

		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		
		nodes.put("X", new Node("X"));
		nodes.put("Y", new Node("Y"));
		
		nodes.put("M", new Node("M"));
		nodes.put("N", new Node("N"));
		
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 0.8));
		
		network.addEdge(new Edge("negative regulation", nodes.get("X"), nodes.get("Y"), false, 0.5));
		
		network.addEdge(new Edge("ptm", nodes.get("M"), nodes.get("N"), false, 0.5, true));

		return network;
	}
	

	/**
	 * Testing the example using console output (use TestFuzzyOverlap for the JUnit version!)
	 * @param args (ignored) argument list
	 */
	public static void main(String[] args)
	{
		FuzzyOverlap ex = new FuzzyOverlap();
		double weight_cutoff = 0.0;
		
		
		System.out.println("Defining network for FuzzyOverlap configuration");
		Project p = ex.getProject();
		int overlap_cutoff = 2;
		//int ID_1 = ex.getTestConfigurationWithReference(p, overlap_cutoff);
		int ID_2 = ex.getTestConfigurationWithoutReference(p, overlap_cutoff);
		
		System.out.print("Calculating 1-all overlap network at weight cutoff " + weight_cutoff);
		System.out.println(" and overlap cutoff " + overlap_cutoff);
		
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID_2, weight_cutoff, -1, 20, true);
		
		System.out.println("");
		ex.printAllNetworks(p, ID_2, false, true);
		
		System.out.println("Log:");
		Logger logger = p.getLogger(ID_2);
		for (LogEntry log : logger.getAllLogMessages())
		{
			System.out.println(log);
		}
	}
	
}
