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
import be.svlandeg.diffany.core.project.LogEntry;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

/**
 * Testing class that simulates input networks with edge conflicts, 
 * e.g. positive and negative regulation at the same time.
 * 
 * @author Sofie Van Landeghem
 */
public class ConflictingEdgesTest extends GenericExample
{

	
	/**
	 * Get a custom project.
	 * @return an example project
	 */
	public Project getDefaultProject()
	{
		String name = "Conflict_test";
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
		ReferenceNetwork r = getTestReference();
		Set<ConditionNetwork> c = getTestCondition();
		boolean cleanInput = true;
		int ID = p.addRunConfiguration(r, c, cleanInput, null);
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
		nodes.put("G", new Node("G", "G"));
		nodes.put("H", new Node("H", "H"));
		nodes.put("J", new Node("J", "J"));
		nodes.put("K", new Node("K", "K"));
		
		ReferenceNetwork network = new ReferenceNetwork("Reference", 1, null);
		
		network.addEdge(new Edge("positive regulation", nodes.get("A"), nodes.get("B"), false, 2, false));
		network.addEdge(new Edge("ptm", nodes.get("A"), nodes.get("B"), false, 5, false));
		network.addEdge(new Edge("somerandomInteraction", nodes.get("A"), nodes.get("B"), false, 4, false));
		
		network.addEdge(new Edge("ptm", nodes.get("G"), nodes.get("H"), false, 3, false));
		network.addEdge(new Edge("phosphorylation", nodes.get("G"), nodes.get("H"), false, 1, false));
		network.addEdge(new Edge("regulation", nodes.get("G"), nodes.get("H"), false, 7, false));
		
		network.addEdge(new Edge("ptm", nodes.get("J"), nodes.get("K"), false, 4, true));
		network.addEdge(new Edge("ubiquitination", nodes.get("J"), nodes.get("K"), true, 2, false));
		network.addEdge(new Edge("regulation", nodes.get("J"), nodes.get("K"), false, 7, false));
		network.addEdge(new Edge("regulation", nodes.get("J"), nodes.get("K"), false, 5, true));
		return network;
	}
	
	/**
	 * Get the custom-defined condition-specific network.
	 * @return the condition-specific network
	 */
	private Set<ConditionNetwork> getTestCondition()
	{
		Set<ConditionNetwork> cnetworks = new HashSet<ConditionNetwork>();

		String description = "Unknown condition";
		Condition c = new Condition(description);
		Set<Condition> conditions = new HashSet<Condition>();
		conditions.add(c);

		ConditionNetwork network = new ConditionNetwork("Condition-specific network", 2, null, conditions);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A", "A"));
		nodes.put("B", new Node("B", "B"));
		nodes.put("G", new Node("G", "G"));
		nodes.put("H", new Node("H", "H"));
		nodes.put("J", new Node("J", "J"));
		nodes.put("K", new Node("K", "K"));
		
		network.addEdge(new Edge("positive regulation", nodes.get("A"), nodes.get("B"), false, 8, false));
		
		network.addEdge(new Edge("phosphorylation", nodes.get("G"), nodes.get("H"), false, 2, false));
		network.addEdge(new Edge("catalysis", nodes.get("G"), nodes.get("H"), false, 4, false));
		
		network.addEdge(new Edge("ptm", nodes.get("J"), nodes.get("K"), true, 3, false));
		network.addEdge(new Edge("ptm", nodes.get("J"), nodes.get("K"), false, 6, false));
		
		network.addEdge(new Edge("positive_regulation", nodes.get("J"), nodes.get("K"), false, 3, true));
		network.addEdge(new Edge("inhibition", nodes.get("J"), nodes.get("K"), true, 4, false));
		
		cnetworks.add(network);
		return cnetworks;
	}
	
	/**
	 * Testing the example
	 * @param args the (ignored) input argument list
	 */
	public static void main(String[] args)
	{
		ConflictingEdgesTest ex = new ConflictingEdgesTest();
		double cutoff = 0.0;
		
		System.out.println("Defining network for conflict test");
		Project p = ex.getDefaultProject();
		int ID = ex.getDefaultRunConfigurationID(p);
		
		System.out.println("Calculating differential networks at cutoff " + cutoff);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, cutoff, true, true, 3, true, null);
		
		System.out.println("");
		ex.printAllNetworks(p, ID, true, false, false);
		
		System.out.println("Logs:");
		for (LogEntry l : p.getLogger(ID).getAllLogMessages())
		{
			System.out.println(l);
		}
	}

}
