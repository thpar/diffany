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
 * This class provides examples taken from the Bandyopadhyay et al, Science 2010 paper.
 * http://www.sciencemag.org/content/330/6009/1385.full.pdf  
 * 
 * @author Sofie Van Landeghem
 */
public class Bandyopadhyay2010 extends GenericExample
{


	/**
	 * Get a custom project.
	 * @return an example project illustrating figure 1C.
	 */
	public Project getDefaultProject()
	{
		String name = "Bandyopadhyay2010_fig1C";
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
		ReferenceNetwork r = getReferenceFigure1C();
		Set<ConditionNetwork> c = getConditionFigure1C();
		boolean cleanInput = true;
		ProgressListener listener = new StandardProgressListener(false);
		int ID = p.addRunConfiguration(r, c, cleanInput, listener);
		return ID;
	}


	/**
	 * Get the reference network depicted in figure 1C.
	 * @return the reference network
	 */
	private ReferenceNetwork getReferenceFigure1C()
	{
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A", "A"));
		nodes.put("B", new Node("B", "B"));
		nodes.put("C", new Node("C", "C"));
		nodes.put("D", new Node("D", "D"));
		nodes.put("E", new Node("E", "E"));
		
		ReferenceNetwork network = new ReferenceNetwork("Untreated Network", 1, null);
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

		ConditionNetwork network = new ConditionNetwork("Treated Network", 2, null, conditions);

		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A", "A"));
		nodes.put("B", new Node("B", "B"));
		nodes.put("C", new Node("C", "C"));
		nodes.put("D", new Node("D", "D"));
		
		network.addEdge(new Edge("negative genetic interaction", nodes.get("A"), nodes.get("D"), true, 1.1));
		network.addEdge(new Edge("positive genetic interaction", nodes.get("A"), nodes.get("B"), true, 0.4));
		network.addEdge(new Edge("negative genetic interaction", nodes.get("A"), nodes.get("C"), true, 0.7));

		cnetworks.add(network);
		return cnetworks;
	}

	/**
	 * Testing the example using console output (use TestExamples for the JUnit version!)
	 * @param args the (ignored) input argument list
	 */
	public static void main(String[] args)
	{
		Bandyopadhyay2010 ex = new Bandyopadhyay2010();
		double cutoff = 0.0;
		
		System.out.println("Defining network for Bandyopadhyay2010 figure 1c");
		Project p = ex.getDefaultProject();
		int ID = ex.getDefaultRunConfigurationID(p);
		ProgressListener listener = new StandardProgressListener(true);
		
		System.out.println("Calculating differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, cutoff, true, true, 3, true, listener);
		
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
