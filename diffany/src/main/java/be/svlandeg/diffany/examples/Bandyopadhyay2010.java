package be.svlandeg.diffany.examples;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.Condition;
import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.Logger;
import be.svlandeg.diffany.concepts.Node;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/** 
 * This class provides examples taken from the Bandyopadhyay et al, Science 2010 paper.
 * http://www.sciencemag.org/content/330/6009/1385.full.pdf  
 * 
 * @author Sofie Van Landeghem
 */
public class Bandyopadhyay2010 extends GenericExample
{

	private NodeMapper nm;
	
	/**
	 * Constructor: generates a default {@link NodeMapper} object
	 */
	public Bandyopadhyay2010()
	{
		nm = new DefaultNodeMapper();
	}
	
	/**
	 * Get a custom project.
	 * @return an example project illustrating figure 1C.
	 */
	public Project getProjectFigure1C()
	{
		String name = "Bandyopadhyay2010_fig1C";
		EdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo, nm);
		return p;
	}
	
	/**
	 * Add some custom-defined networks to the project.
	 * @return the resulting configuration ID.
	 */
	public int getTestConfiguration1C(Project p)
	{
		ReferenceNetwork r = getReferenceFigure1C();
		Set<ConditionNetwork> c = getConditionFigure1C();
		int ID = p.addRunConfiguration(r, c);
		return ID;
	}


	/**
	 * Get the reference network depicted in figure 1C.
	 * @return the reference network
	 */
	private ReferenceNetwork getReferenceFigure1C()
	{
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("C", new Node("C"));
		nodes.put("D", new Node("D"));
		nodes.put("E", new Node("E"));
		
		ReferenceNetwork network = new ReferenceNetwork("Untreated Network", nm);
		network.addEdge(new Edge("negative genetic interaction", nodes.get("A"), nodes.get("D"), true, 1.1));
		network.addEdge(new Edge("negative genetic interaction", nodes.get("A"), nodes.get("B"), true, 0.3));
		network.addEdge(new Edge("positive genetic interaction", nodes.get("E"), nodes.get("C"), true, 0.8));
		return network;
	}

	/**
	 * Get the condition-specific network depicted in figure 1C.
	 * @return the condition-specific network
	 */
	private Set<ConditionNetwork> getConditionFigure1C()
	{
		Set<ConditionNetwork> cnetworks = new HashSet<ConditionNetwork>();

		String description = "treated with MMS";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Treated Network", conditions, nm);

		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("C", new Node("C"));
		nodes.put("D", new Node("D"));
		
		network.addEdge(new Edge("negative genetic interaction", nodes.get("A"), nodes.get("D"), true, 1.1));
		network.addEdge(new Edge("positive genetic interaction", nodes.get("A"), nodes.get("B"), true, 0.4));
		network.addEdge(new Edge("negative genetic interaction", nodes.get("A"), nodes.get("C"), true, 0.7));

		cnetworks.add(network);
		return cnetworks;
	}

	/**
	 * Testing the example using console output (use TestExamples for the JUnit version!)
	 */
	public static void main(String[] args)
	{
		Bandyopadhyay2010 ex = new Bandyopadhyay2010();
		double cutoff = 0.0;
		
		System.out.println("Defining network for Bandyopadhyay2010 figure 1c");
		Project p = ex.getProjectFigure1C();
		int ID = ex.getTestConfiguration1C(p);
		
		System.out.println("Calculating differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, cutoff);
		
		System.out.println("");
		ex.printAllNetworks(p, ID);
		
		System.out.println("Log:");
		Logger logger = p.getLogger(ID);
		for (String log : logger.getAllLogMessages())
		{
			System.out.println(log);
		}
	}
	
}
