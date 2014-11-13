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
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

/** 
 * This class provides examples to benchmark the fuzzy consensus & differential functionality.  
 * Specifically, we use the same interaction type for all edges, but with different weights and complex differences between various cutoffs.
 * 
 * @author Sofie Van Landeghem
 */
public class FuzzyNetworks2 extends GenericExample
{

	/**
	 * Constructor
	 */
	public FuzzyNetworks2()
	{}
	
	/**
	 * Get a custom project.
	 * @return an example project illustrating figure 1C.
	 */
	public Project getDefaultProject()
	{
		String name = "FuzzyNetworks2";
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo);
		return p;
	}
	
	/**
	 * Add some custom-defined networks to the project: 1 reference network and 3 condition-specific.
	 * @param p the fuzzy project
	 * @return the resulting configuration ID.
	 */
	public int getDefaultRunConfigurationID(Project p)
	{
		return getDefaultRunConfigurationID(p, 4);
	}
	
	/**
	 * Add some custom-defined networks to the project: 1 reference network and 3 condition-specific.
	 * @param p the fuzzy project
	 * @param supportingCutoff the number of networks that should at least match for consensus to be defined: min. 2, max conditions.size + 1.
	 * @return the resulting configuration ID.
	 */
	public int getDefaultRunConfigurationID(Project p, int supportingCutoff)
	{
		ReferenceNetwork r = getReference();
		Set<ConditionNetwork> c = new HashSet<ConditionNetwork>();
		c.add(getCondition1());
		c.add(getCondition2());
		c.add(getCondition3());	
		boolean cleanInput = true;
		int ID = p.addRunConfiguration(r, c, supportingCutoff, cleanInput, null);
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
		
		ReferenceNetwork network = new ReferenceNetwork("Fuzzy2 reference network", 1, null);
		
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 0.8, false));
		network.addEdge(new Edge("regulates", nodes.get("X"), nodes.get("Y"), false, 11, false));
		
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

		ConditionNetwork network = new ConditionNetwork("Condition network 1", 11, null, conditions);

		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		
		nodes.put("X", new Node("X"));
		nodes.put("Y", new Node("Y"));
		
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 1.2, false));

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

		ConditionNetwork network = new ConditionNetwork("Condition network 2", 12, null, conditions);

		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		
		nodes.put("X", new Node("X"));
		nodes.put("Y", new Node("Y"));
		
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 0.6, false));
		network.addEdge(new Edge("regulates", nodes.get("X"), nodes.get("Y"), false, 6, false));

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

		ConditionNetwork network = new ConditionNetwork("Condition network 3", 13, null, conditions);

		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		
		nodes.put("X", new Node("X"));
		nodes.put("Y", new Node("Y"));
		
		network.addEdge(new Edge("ppi", nodes.get("B"), nodes.get("A"), true, 0.4, false));
		network.addEdge(new Edge("regulates", nodes.get("X"), nodes.get("Y"), false, 4, false));

		return network;
	}
	
	
	/**
	 * Testing the example using console output (use TestFuzzyConsensus for the JUnit version!)
	 * @param args (ignored) argument list
	 */
	public static void main(String[] args)
	{
		FuzzyNetworks2 ex = new FuzzyNetworks2();
		double weightCutoff = 0.0;
		
		System.out.println("Defining network for FuzzyNetworks2 configuration");
		Project p = ex.getDefaultProject();
		int supportingCutoff = 4;
		int ID = ex.getDefaultRunConfigurationID(p, supportingCutoff);
		
		System.out.print("Calculating 1-all differential networks at weight cutoff " + weightCutoff);
		System.out.println(" and supporting cutoff " + supportingCutoff);
		
		boolean minOperator = true;
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weightCutoff, 70, 80, minOperator, null);	
		
		System.out.println("");
		ex.printAllNetworks(p, ID, true, false, false);	
		
		System.out.println("Log:");
		Logger logger = p.getLogger(ID);
		for (LogEntry log : logger.getAllLogMessages())
		{
			System.out.println(log);
		}
	}
	
}
