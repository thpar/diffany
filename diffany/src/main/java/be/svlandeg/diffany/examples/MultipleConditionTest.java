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
 * Testing class that tries to simulate a use-case of calculating a differential network
 * from one reference network and multiple condition-specific networks.
 * 
 * @author Sofie Van Landeghem
 */
public class MultipleConditionTest extends GenericExample
{
	
	private NodeMapper nm;
	
	/**
	 * Constructor: generates a default {@link NodeMapper} object
	 */
	public MultipleConditionTest()
	{
		nm = new DefaultNodeMapper();
	}
	
	/**
	 * Get a custom project.
	 * @return an example project
	 */
	public Project getTestProject()
	{
		String name = "Multiple_test";
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo, nm);
		return p;
	}
	
	/**
	 * Add some custom-defined networks to the project.
	 * @param p the input project
	 * @return the resulting configuration ID
	 */
	public int getTestDiffConfiguration(Project p)
	{
		ReferenceNetwork r = getTestReference();
		Set<ConditionNetwork> c = getTestConditions();
		int ID = p.addRunConfiguration(r, c);
		return ID;
	}
	
	/**
	 * Add some custom-defined networks to the project.
	 * @param p the input project
	 * @return the resulting configuration ID
	 */
	public int getTestConsensusConfiguration(Project p)
	{
		Set<InputNetwork> i = new HashSet<InputNetwork>();
		i.addAll(getTestConditions());
		i.add(getTestReference());
		int ID = p.addRunConfiguration(i);
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
		nodes.put("C", new Node("C"));
		nodes.put("W", new Node("W"));
		nodes.put("Z", new Node("Z"));
		
		ReferenceNetwork network = new ReferenceNetwork("Reference", 1, nm);
		
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 0.7, false));
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("C"), true, 0.8, false));
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("Z"), true, 0.9, false));
		network.addEdge(new Edge("ppi", nodes.get("Z"), nodes.get("W"), true, 0.5, false));
		
		network.addEdge(new Edge("phosphorylates", nodes.get("A"), nodes.get("B"), false, 1.1, false));
		
		nodes.put("M", new Node("M"));
		nodes.put("N", new Node("N"));
		nodes.put("O", new Node("O"));
		nodes.put("P", new Node("P"));
		network.addEdge(new Edge("phosphorylates", nodes.get("M"), nodes.get("N"), false, 2, false));
		network.addEdge(new Edge("phosphorylation", nodes.get("M"), nodes.get("O"), false, 8, true));
		network.addEdge(new Edge("phosphorylation", nodes.get("M"), nodes.get("P"), false, 4, true));
		network.addEdge(new Edge("phosphorylation", nodes.get("N"), nodes.get("P"), false, 3, false));
		network.addEdge(new Edge("phosphorylation", nodes.get("N"), nodes.get("O"), false, 1, true));
		network.addEdge(new Edge("ptm", nodes.get("O"), nodes.get("P"), false, 8, true));

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

		ConditionNetwork network = new ConditionNetwork("Salty", 2, conditions, nm);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("C", new Node("C"));
		nodes.put("D", new Node("D"));
		nodes.put("F", new Node("F"));
		nodes.put("W", new Node("W"));
		nodes.put("Z", new Node("Z"));

		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 0.4, false));
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("C"), true, 0.6, false));
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("Z"), true, 0.1, false));
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("D"), true, 0.9, false));
		network.addEdge(new Edge("ppi", nodes.get("D"), nodes.get("F"), true, 0.3, false));
		
		nodes.put("M", new Node("M"));
		nodes.put("N", new Node("N"));
		nodes.put("O", new Node("O"));
		nodes.put("P", new Node("P"));
		network.addEdge(new Edge("phosphorylates", nodes.get("M"), nodes.get("N"), false, 8, false));
		network.addEdge(new Edge("phosphorylation", nodes.get("O"), nodes.get("N"), false, 5, true));
		network.addEdge(new Edge("phosphorylates", nodes.get("O"), nodes.get("P"), false, 4, false));
		network.addEdge(new Edge("phosphorylation", nodes.get("P"), nodes.get("M"), false, 2, false));
		network.addEdge(new Edge("phosphorylation", nodes.get("M"), nodes.get("O"), false, 1, true));

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

		ConditionNetwork network = new ConditionNetwork("Draughty", 3, conditions, nm);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("C", new Node("C"));
		nodes.put("D", new Node("D"));
		nodes.put("E", new Node("E"));
		nodes.put("W", new Node("W"));
		nodes.put("Z", new Node("Z"));

		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 0.3, false));
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("C"), true, 1.2, false));
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("D"), true, 0.75, false));
		network.addEdge(new Edge("ppi", nodes.get("D"), nodes.get("E"), true, 0.2, false));
				
		nodes.put("M", new Node("M"));
		nodes.put("N", new Node("N"));
		nodes.put("O", new Node("O"));
		nodes.put("P", new Node("P"));
		network.addEdge(new Edge("phosphorylate", nodes.get("M"), nodes.get("N"), false, 6, false));
		network.addEdge(new Edge("phosphorylation", nodes.get("N"), nodes.get("P"), false, 4, true));
		network.addEdge(new Edge("phosphorylation", nodes.get("N"), nodes.get("O"), false, 5, false));
		network.addEdge(new Edge("ubiquitinates", nodes.get("P"), nodes.get("O"), false, 6, true));
		network.addEdge(new Edge("ptm", nodes.get("P"), nodes.get("M"), false, 7, false));
		network.addEdge(new Edge("phosphorylation", nodes.get("P"), nodes.get("N"), false, 8, false));

		return network;
	}
	
	/**
	 * Testing the example
	 * @param args the (ignored) input argument list
	 */
	public static void main(String[] args)
	{
		MultipleConditionTest ex = new MultipleConditionTest();
		double cutoff = 0.0;
		
		System.out.println("Defining network for weight test");
		Project p = ex.getTestProject();
		int ID_diff = ex.getTestDiffConfiguration(p);
		
		System.out.println("Calculating 1-all differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID_diff, cutoff, 10, 11, true);
		
		System.out.println("");
		ex.printAllNetworks(p, ID_diff, true, false, false);
		
		System.out.println(" **************************************************************** ");
		
		Logger log = p.getLogger(ID_diff);
		for (LogEntry e: log.getAllLogMessages())
		{
			System.out.println(e);
		}
		
		/*
		System.out.println(" ");
		System.out.println(" **************************************************************** ");
		System.out.println(" ");
		
		System.out.println("Calculating pairwise differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID_diff, cutoff, true, true);
		
		System.out.println("");
		ex.printAllNetworks(p, ID_diff);
		
		System.out.println(" ");
		System.out.println(" **************************************************************** ");
		System.out.println(" ");
		
		System.out.println("Calculating pairwise consensus networks at cutoff " + cutoff);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID_diff, cutoff, false, true);
		
		System.out.println("");
		ex.printAllConsensusNetworks(p, ID_diff);
		*/
		
	}
}
