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
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

/**
 * Testing class that simulates a few of the osmotic use-case edges (4 time-specific condition networks).
 * 
 * @author Sofie Van Landeghem
 */
public class OsmoticSampleTest extends GenericExample
{
	
	
	/**
	 * Get a custom project.
	 * @return an example project
	 */
	public Project getDefaultProject()
	{
		String name = "Osmotic_sample_test";
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
		return getTestDiffConfiguration(p, 5);
	}
	
	/**
	 * Add some custom-defined networks to the project.
	 * @param p the input project
	 * @param supportingCutoff the number of networks that should at least match for consensus to be defined: min. 2, max conditions.size + 1.
	 * @return the resulting configuration ID
	 */
	public int getTestDiffConfiguration(Project p, int supportingCutoff)
	{
		ReferenceNetwork r = getTestReference();
		Set<ConditionNetwork> c = getTestConditions();
		boolean cleanInput = true;
		int ID = p.addRunConfiguration(r, c, supportingCutoff, cleanInput, null);
		return ID;
	}
	
	/**
	 * Add some custom-defined networks to the project.
	 * @param p the input project
	 * @param supportingCutoff the number of networks that should at least match for consensus to be defined.
	 * @return the resulting configuration ID
	 */
	public int getTestConsensusConfiguration(Project p, int supportingCutoff)
	{
		Set<InputNetwork> i = new HashSet<InputNetwork>();
		i.addAll(getTestConditions());
		i.add(getTestReference());
		int ID = p.addRunConfiguration(i, supportingCutoff, false, null);
		return ID;
	}

	/**
	 * Get a custom-defined reference network.
	 * @return the reference network
	 */
	private ReferenceNetwork getTestReference()
	{
		Map<String, Node> nodes = new HashMap<String, Node>();
		
		nodes.put("A", new Node("A", "A"));
		nodes.put("B", new Node("B", "B"));
		nodes.put("X", new Node("X", "X"));
		nodes.put("Y", new Node("Y", "Y"));
		
		nodes.put("M", new Node("M", "M"));
		nodes.put("N", new Node("N", "N"));
		nodes.put("O", new Node("O", "O"));
		nodes.put("P", new Node("P", "P"));
		
		ReferenceNetwork network = new ReferenceNetwork("Reference", 1, null);
		
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 1.0, false));
		network.addEdge(new Edge("ppi", nodes.get("X"), nodes.get("Y"), true, 1.0, false));
		
		network.addEdge(new Edge("downregulation", nodes.get("M"), nodes.get("N"), false, 0, false));
		network.addEdge(new Edge("downregulation", nodes.get("O"), nodes.get("P"), false, 1.0, false));

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

		ConditionNetwork network = new ConditionNetwork("Mannitol after 1.5h", 2, null, conditions);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		
		nodes.put("A", new Node("A", "A"));
		nodes.put("B", new Node("B", "B"));
		nodes.put("X", new Node("X", "X"));
		nodes.put("Y", new Node("Y", "Y"));
		
		nodes.put("M", new Node("M", "M"));
		nodes.put("N", new Node("N", "N"));
		nodes.put("O", new Node("O", "O"));
		nodes.put("P", new Node("P", "P"));
		
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 1.0, false));
		network.addEdge(new Edge("ppi", nodes.get("X"), nodes.get("Y"), true, 1.0, false));
		
		network.addEdge(new Edge("downregulation", nodes.get("M"), nodes.get("N"), false, 1.2, false));
		network.addEdge(new Edge("downregulation", nodes.get("O"), nodes.get("P"), false, 2.2, false));

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

		ConditionNetwork network = new ConditionNetwork("Mannitol after 3h", 3, null, conditions);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		
		nodes.put("A", new Node("A", "A"));
		nodes.put("B", new Node("B", "B"));
		nodes.put("X", new Node("X", "X"));
		nodes.put("Y", new Node("Y", "Y"));
		
		nodes.put("M", new Node("M", "M"));
		nodes.put("N", new Node("N", "N"));
		nodes.put("O", new Node("O", "O"));
		nodes.put("P", new Node("P", "P"));
		
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 1.3, false));
		network.addEdge(new Edge("ppi", nodes.get("X"), nodes.get("Y"), true, 1.3, false));
		
		network.addEdge(new Edge("downregulation", nodes.get("M"), nodes.get("N"), false, 1.4, false));
		network.addEdge(new Edge("downregulation", nodes.get("O"), nodes.get("P"), false, 2.4, false));

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

		ConditionNetwork network = new ConditionNetwork("Mannitol after 12h", 4, null, conditions);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		
		nodes.put("A", new Node("A", "A"));
		nodes.put("B", new Node("B", "B"));
		nodes.put("X", new Node("X", "X"));
		nodes.put("Y", new Node("Y", "Y"));
		
		nodes.put("M", new Node("M", "M"));
		nodes.put("N", new Node("N", "N"));
		nodes.put("O", new Node("O", "O"));
		nodes.put("P", new Node("P", "P"));
		
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 1.3, false));
		network.addEdge(new Edge("ppi", nodes.get("X"), nodes.get("Y"), true, 1.3, false));
		
		network.addEdge(new Edge("downregulation", nodes.get("M"), nodes.get("N"), false, 1.5, false));
		network.addEdge(new Edge("downregulation", nodes.get("O"), nodes.get("P"), false, 2.5, false));

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

		ConditionNetwork network = new ConditionNetwork("Mannitol after 24h", 5, null, conditions);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		
		nodes.put("A", new Node("A", "A"));
		nodes.put("B", new Node("B", "B"));
		nodes.put("X", new Node("X", "X"));
		nodes.put("Y", new Node("Y", "Y"));
		
		nodes.put("M", new Node("M", "M"));
		nodes.put("N", new Node("N", "N"));
		nodes.put("O", new Node("O", "O"));
		nodes.put("P", new Node("P", "P"));
		
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 0.0, false));
		// no X-Y edge
		
		network.addEdge(new Edge("downregulation", nodes.get("M"), nodes.get("N"), false, 1.1, false));
		network.addEdge(new Edge("downregulation", nodes.get("O"), nodes.get("P"), false, 2.1, false));

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
		Project p = ex.getDefaultProject();
		int ID_diff = ex.getTestDiffConfiguration(p, 5);
		
		System.out.println("Calculating 1-all differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID_diff, cutoff, null, null, 10, 11, true, null);
		
		System.out.println("");
		ex.printAllNetworks(p, ID_diff, true, false, false);
		
		System.out.println(" **************************************************************** ");
		
		Logger log = p.getLogger(ID_diff);
		for (LogEntry e: log.getAllLogMessages())
		{
			System.out.println(e);
		}
	}

}
