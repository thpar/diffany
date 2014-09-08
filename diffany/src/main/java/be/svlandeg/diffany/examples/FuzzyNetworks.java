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
 * This class provides examples to benchmark the fuzzy consensus functionality.  
 * Specifically, different interaction types, weights and negation states are mixed.
 * 
 * @author Sofie Van Landeghem
 */
public class FuzzyNetworks extends GenericExample
{

	private NodeMapper nm;
	
	/**
	 * Constructor: generates a default {@link NodeMapper} object
	 */
	public FuzzyNetworks()
	{
		nm = new DefaultNodeMapper();
	}
	
	/**
	 * Get a custom project.
	 * @return an example project illustrating figure 1C.
	 */
	public Project getProject()
	{
		String name = "FuzzyNetworks";
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo, nm);
		return p;
	}
	
	/**
	 * Add some custom-defined networks to the project: 1 reference network and 3 condition-specific.
	 * @param p the fuzzy project
	 * @param supportingCutoff the minimal number of networks that need to agree on a certain edge
	 * @return the resulting configuration ID.
	 */
	public int getTestConfigurationWithReference(Project p, int supportingCutoff)
	{
		ReferenceNetwork r = getReference();
		Set<ConditionNetwork> c = new HashSet<ConditionNetwork>();
		c.add(getCondition1());
		c.add(getCondition2());
		c.add(getCondition3());	
		int ID = p.addRunConfiguration(r, c, supportingCutoff);
		return ID;
	}
	
	/**
	 * Add some custom-defined networks to the project: 4 condition-specific networks, no reference.
	 * 
	 * @param p the fuzzy project
	 * @param supportingCutoff the minimal number of networks that need to agree on a certain edge
	 * @param refAsCond whether or not to include the reference network as "condition 0". If not, only 3 input networks are used.
	 * @return the resulting configuration ID.
	 */
	public int getTestConfigurationWithoutReference(Project p, int supportingCutoff, boolean refAsCond)
	{
		Set<InputNetwork> c = new HashSet<InputNetwork>();
		if (refAsCond)
		{
			c.add(getCondition0());	
		}
		c.add(getCondition1());	
		c.add(getCondition2());
		c.add(getCondition3());	
		int ID = p.addRunConfiguration(c, supportingCutoff, false);
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
		
		network.addEdge(new Edge("colocalization", nodes.get("A"), nodes.get("B"), false, 0.6, false));
		
		network.addEdge(new Edge("negative regulation", nodes.get("X"), nodes.get("Y"), false, 0.5, false));
		
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
		
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 0.8));
		
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
		
		network.addEdge(new Edge("negative regulation", nodes.get("X"), nodes.get("Y"), false, 0.6));
		
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
		
		network.addEdge(new Edge("positive regulation", nodes.get("X"), nodes.get("Y"), false, 0.7));
		
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
		
		network.addEdge(new Edge("colocalization", nodes.get("A"), nodes.get("B"), false, 0.6));
		
		network.addEdge(new Edge("negative regulation", nodes.get("X"), nodes.get("Y"), false, 0.5));
		
		network.addEdge(new Edge("ptm", nodes.get("M"), nodes.get("N"), false, 0.5, true));

		return network;
	}
	

	/**
	 * Testing the example using console output (use TestFuzzyConsensus for the JUnit version!)
	 * @param args (ignored) argument list
	 */
	public static void main(String[] args)
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		
		System.out.println("Defining network for FuzzyNetworks configuration");
		Project p = ex.getProject();
		int supportingCutoff = 2;
		int ID = ex.getTestConfigurationWithReference(p, supportingCutoff);
		//int ID = ex.getTestConfigurationWithoutReference(p, supportingCutoff, true);
		
		System.out.print("Calculating 1-all consensus network at weight cutoff " + weight_cutoff);
		System.out.println(" and supporting networks cutoff " + supportingCutoff);
		
		//new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 20, true);	// consensus, no diff
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, 20, -1, true, null);	// diff, no consensus
		
		System.out.println("");
		//ex.printAllNetworks(p, ID, false, true, false);	// consensus, no diff
		ex.printAllNetworks(p, ID, false, false, true);	// diff, no consensus
		
		System.out.println("Log:");
		Logger logger = p.getLogger(ID);
		for (LogEntry log : logger.getAllLogMessages())
		{
			System.out.println(log);
		}
	}
	
}
