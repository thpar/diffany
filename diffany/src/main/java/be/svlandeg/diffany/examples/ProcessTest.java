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

import java.util.*;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.*;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.semantics.*;

/**
 * Testing class that tries to simulate a range of possibilities in process networks
 * and is able to produce differential networks for them.
 * 
 * @author Sofie Van Landeghem
 */
public class ProcessTest extends GenericExample
{
	
	
	/**
	 * Get a custom project.
	 * @return an example project
	 */
	public Project getDefaultProject()
	{
		String name = "Process_test";
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
		nodes.put("M", new Node("M", "M"));
		nodes.put("N", new Node("N", "N"));
		nodes.put("S", new Node("S", "S"));
		nodes.put("T", new Node("T", "T"));
		nodes.put("X", new Node("X", "X"));
		nodes.put("Y", new Node("Y", "Y"));
		nodes.put("G", new Node("G", "G"));
		nodes.put("H", new Node("H", "H"));
		nodes.put("J", new Node("J", "J"));
		nodes.put("K", new Node("K", "K"));
		
		ReferenceNetwork network = new ReferenceNetwork("Condition 1", 1, null);
		network.addEdge(new Edge("ppi", nodes.get("A"), nodes.get("B"), true, 2, false));
		network.addEdge(new Edge("ppi", nodes.get("S"), nodes.get("T"), true, 3, false));
		network.addEdge(new Edge("ptm", nodes.get("X"), nodes.get("Y"), true, 4, false));
		network.addEdge(new Edge("ubiquitinates", nodes.get("H"), nodes.get("G"), true, 1, false));
		network.addEdge(new Edge("ptm", nodes.get("J"), nodes.get("K"), true, 5, true));
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

		ConditionNetwork network = new ConditionNetwork("Condition 2", 2, null, conditions);
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		nodes.put("A", new Node("A", "A"));
		nodes.put("B", new Node("B", "B"));
		nodes.put("M", new Node("M", "M"));
		nodes.put("N", new Node("N", "N"));
		nodes.put("S", new Node("S", "S"));
		nodes.put("T", new Node("T", "T"));
		nodes.put("X", new Node("X", "X"));
		nodes.put("Y", new Node("Y", "Y"));
		nodes.put("G", new Node("G", "G"));
		nodes.put("H", new Node("H", "H"));
		nodes.put("J", new Node("J", "J"));
		nodes.put("K", new Node("K", "K"));

		network.addEdge(new Edge("ppi", nodes.get("M"), nodes.get("N"), true, 3, false));
		network.addEdge(new Edge("phosphorylates", nodes.get("S"), nodes.get("T"), true, 2, false));
		network.addEdge(new Edge("phosphorylates", nodes.get("X"), nodes.get("Y"), true, 3, false));
		network.addEdge(new Edge("phosphorylation", nodes.get("G"), nodes.get("H"), false, 5, false));
		network.addEdge(new Edge("phosphorylate", nodes.get("K"), nodes.get("J"), true, 4, true));
		
		cnetworks.add(network);
		return cnetworks;
	}
	
	/**
	 * Testing the example
	 * @param args the (ignored) input argument list
	 */
	public static void main(String[] args)
	{
		ProcessTest ex = new ProcessTest();
		double cutoff = 0.0;
		
		System.out.println("Defining network for weight test");
		Project p = ex.getDefaultProject();
		int ID = ex.getDefaultRunConfigurationID(p);
		
		System.out.println("Calculating differential networks at cutoff " + cutoff);
		//new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, cutoff, true, true, 3, true, null);
		//ex.printAllNetworks(p, ID, true, false, false);
		
		/* Test for bug issue #233 */
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, cutoff, null, null, 5, 3, true, null);
		ex.printAllNetworks(p, ID, false, true, false);
		
		System.out.println("");
		
	}

}
