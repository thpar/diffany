package be.svlandeg.diffany.examples;

/*
 * #%L
 * Diffany
 * %%
 * Copyright (C) 2014 PSB/UGent - Sofie Van Landeghem and Thomas Van Parys
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU General Lesser Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/lgpl-3.0.html>.
 * #L%
 */

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
import be.svlandeg.diffany.core.progress.ProgressListener;
import be.svlandeg.diffany.core.progress.StandardProgressListener;
import be.svlandeg.diffany.core.project.LogEntry;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
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
	
	/**
	 * Constructor
	 */
	public Ideker2011()
	{}
	
	/**
	 * Get a custom project.
	 * @return an example project illustrating figure 1C
	 */
	public Project getDefaultProject()
	{
		String name = "Ideker2011_fig3A";
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo);
		return p;
	}
	
	/**
	 * Add some custom-defined networks to the project.
	 * @param p the input project
	 * @return the resulting configuration ID
	 */
	public int getDefaultRunConfigurationID(Project p)
	{
		ReferenceNetwork r = getReferenceFigure3A();
		ProgressListener listener = new StandardProgressListener(false);
		
		Set<ConditionNetwork> c = getConditionFigure3A();
		boolean cleanInput = false;
		int ID = p.addRunConfiguration(r, c, cleanInput, listener);
		
		return ID;
	}

	/**
	 * Get the reference network depicted in figure 3A.
	 * @return the reference network
	 */
	private ReferenceNetwork getReferenceFigure3A()
	{
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A", "A"));
		nodes.put("B", new Node("B", "B"));
		nodes.put("D", new Node("D", "D"));
		nodes.put("E", new Node("E", "E"));
		nodes.put("F", new Node("F", "F"));
		
		ReferenceNetwork network = new ReferenceNetwork("Reference Network Ideker2011", 1, null);
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

		ConditionNetwork network = new ConditionNetwork("Condition Network Ideker2011", 2, null, conditions);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A", "A"));
		nodes.put("B", new Node("B", "B"));
		nodes.put("C", new Node("C", "C"));
		nodes.put("D", new Node("D", "D"));
		nodes.put("F", new Node("F", "F"));

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
		Project p = ex.getDefaultProject();
		int ID = ex.getDefaultRunConfigurationID(p);
		
		ProgressListener listener = new StandardProgressListener(true);
		
		System.out.println("Calculating differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, cutoff, true, true, 3, true, listener);

		System.out.println("");
		ex.printAllNetworks(p, ID, true, false, false);
		
		Logger l = p.getLogger(ID);
		for (LogEntry msg : l.getAllLogMessages())
		{
			System.out.println(msg);
		}
	}
}
