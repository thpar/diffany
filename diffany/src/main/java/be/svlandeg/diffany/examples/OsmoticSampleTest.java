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
 * Testing class that simulates a few of the osmotic use-case edges (4 time-specific condition networks).
 * 
 * @author Sofie Van Landeghem
 */
public class OsmoticSampleTest extends GenericExample
{
private NodeMapper nm;
	
	/**
	 * Constructor: generates a default {@link NodeMapper} object
	 */
	public OsmoticSampleTest()
	{
		nm = new DefaultNodeMapper();
	}
	
	/**
	 * Get a custom project.
	 * @return an example project
	 */
	public Project getTestProject()
	{
		String name = "Osmotic_sample_test";
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
		boolean cleanInput = true;
		int ID = p.addRunConfiguration(r, c, cleanInput);
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
		
		nodes.put("at1g06225", new Node("at1g06225"));
		nodes.put("at1g65690", new Node("at1g65690"));
		
		ReferenceNetwork network = new ReferenceNetwork("Reference", 1, nm);
		
		network.addEdge(new Edge("ppi", nodes.get("at1g06225"), nodes.get("at1g65690"), true, 1.0, false));

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
		cnetworks.add(getThirdCondition());
		cnetworks.add(getFourthCondition());
		return cnetworks;
	}
	
	/**
	 * Get the first custom-defined condition-specific network.
	 * @return the first condition-specific network
	 */
	private ConditionNetwork getFirstCondition()
	{
		String description = "1.5h";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Mannitol after 1.5h", 2, conditions, nm);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		
		nodes.put("at1g06225", new Node("at1g06225"));
		nodes.put("at1g65690", new Node("at1g65690"));
		
		network.addEdge(new Edge("ppi", nodes.get("at1g06225"), nodes.get("at1g65690"), true, 1.0, false));

		return network;
	}
	
	/**
	 * Get the second custom-defined condition-specific network.
	 * @return the second condition-specific network
	 */
	private ConditionNetwork getSecondCondition()
	{
		String description = "3h";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Mannitol after 3h", 3, conditions, nm);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		
		nodes.put("at1g06225", new Node("at1g06225"));
		nodes.put("at1g65690", new Node("at1g65690"));
		
		network.addEdge(new Edge("ppi", nodes.get("at1g06225"), nodes.get("at1g65690"), true, 1.3, false));

		return network;
	}
	
	/**
	 * Get the third custom-defined condition-specific network.
	 * @return the third condition-specific network
	 */
	private ConditionNetwork getThirdCondition()
	{
		String description = "12h";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Mannitol after 12h", 4, conditions, nm);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		
		nodes.put("at1g06225", new Node("at1g06225"));
		nodes.put("at1g65690", new Node("at1g65690"));
		
		network.addEdge(new Edge("ppi", nodes.get("at1g06225"), nodes.get("at1g65690"), true, 1.3, false));

		return network;
	}
	
	/**
	 * Get the fourth custom-defined condition-specific network.
	 * @return the fourth condition-specific network
	 */
	private ConditionNetwork getFourthCondition()
	{
		String description = "24h";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Mannitol after 24h", 5, conditions, nm);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		
		nodes.put("at1g06225", new Node("at1g06225"));
		nodes.put("at1g65690", new Node("at1g65690"));
		
		network.addEdge(new Edge("ppi", nodes.get("at1g06225"), nodes.get("at1g65690"), true, 0.0, false));

		return network;
	}
	
	/**
	 * Testing the example
	 * @param args the (ignored) input argument list
	 */
	public static void main(String[] args)
	{
		OsmoticSampleTest ex = new OsmoticSampleTest();
		double cutoff = 0.0;
		
		System.out.println("Defining network for weight test");
		Project p = ex.getTestProject();
		int ID_diff = ex.getTestDiffConfiguration(p);
		
		System.out.println("Calculating 1-all differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID_diff, cutoff, 10, 11, true, null);
		
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
