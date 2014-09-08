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
import be.svlandeg.diffany.core.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;


/** 
 * This class provides examples taken from the Ideker et al, Molecular Systems Biology 2011 paper.
 * http://www.nature.com/msb/journal/v8/n1/pdf/msb201199.pdf
 * 
 * In contrast to the original picture, the Diffany algorithm will also create a neutral
 * 'regulation' edge when positive regulation turns into negative regulation or vice versa.
 * 
 * @author Sofie Van Landeghem
 */
public class Ideker2011 extends GenericExample
{

	private NodeMapper nm;
	
	/**
	 * Constructor: generates a default {@link NodeMapper} object
	 */
	public Ideker2011()
	{
		nm = new DefaultNodeMapper();
	}
	
	/**
	 * Get a custom project.
	 * @return an example project illustrating figure 1C
	 */
	public Project getProjectFigure3A()
	{
		String name = "Ideker2011_fig3A";
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo, nm);
		return p;
	}
	
	/**
	 * Add some custom-defined networks to the project.
	 * @param p the input project
	 * @return the resulting configuration ID
	 */
	public int getTestConfiguration3A(Project p)
	{
		ReferenceNetwork r = getReferenceFigure3A();
		Set<ConditionNetwork> c = getConditionFigure3A();
		int ID = p.addRunConfiguration(r, c);
		return ID;
	}

	/**
	 * Get the reference network depicted in figure 3A.
	 * @return the reference network
	 */
	private ReferenceNetwork getReferenceFigure3A()
	{
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("D", new Node("D"));
		nodes.put("E", new Node("E"));
		nodes.put("F", new Node("F"));
		
		ReferenceNetwork network = new ReferenceNetwork("Reference Network Ideker2011", 1, nm);
		network.addEdge(new Edge("negative genetic interaction", nodes.get("A"), nodes.get("F"), true, 1.0));
		network.addEdge(new Edge("negative genetic interaction", nodes.get("A"), nodes.get("D"), true, 0.7));
		network.addEdge(new Edge("negative genetic interaction", nodes.get("A"), nodes.get("B"), true, 0.3));
		network.addEdge(new Edge("positive genetic interaction", nodes.get("A"), nodes.get("E"), true, 0.8));
		return network;
	}

	/**
	 * Get the condition-specific network depicted in figure 3A.
	 * @return the condition-specific network
	 */
	private Set<ConditionNetwork> getConditionFigure3A()
	{
		Set<ConditionNetwork> cnetworks = new HashSet<ConditionNetwork>();

		String description = "Unknown condition";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Condition Network Ideker2011", 2, conditions, nm);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A"));
		nodes.put("B", new Node("B"));
		nodes.put("D", new Node("D"));
		nodes.put("C", new Node("C"));
		nodes.put("F", new Node("F"));

		network.addEdge(new Edge("negative genetic interaction", nodes.get("F"), nodes.get("A"), true, 1.0));
		network.addEdge(new Edge("positive genetic interaction", nodes.get("B"), nodes.get("A"), true, 0.4));
		network.addEdge(new Edge("negative genetic interaction", nodes.get("C"), nodes.get("A"), true, 1.2));
		network.addEdge(new Edge("negative genetic interaction", nodes.get("D"), nodes.get("A"), true, 0.7));
		
		cnetworks.add(network);
		return cnetworks;
	}

	/**
	 * Testing the example using console output (use TestExamples for the JUnit version!)
	 * @param args the (ignored) input argument list
	 */
	public static void main(String[] args)
	{
		Ideker2011 ex = new Ideker2011();
		double cutoff = 0.0;
		
		System.out.println("Defining network for Ideker2011 figure 3A");
		Project p = ex.getProjectFigure3A();
		int ID = ex.getTestConfiguration3A(p);
		
		System.out.println("Calculating differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, cutoff, true, true, 3, true, null);
		
		System.out.println("");
		ex.printAllNetworks(p, ID, true, false, false);
		
		Logger l = p.getLogger(ID);
		for (LogEntry msg : l.getAllLogMessages())
		{
			System.out.println(msg);
		}
	}
}
