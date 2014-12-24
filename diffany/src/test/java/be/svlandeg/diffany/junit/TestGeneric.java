package be.svlandeg.diffany.junit;

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


import static org.junit.Assert.assertEquals;

import java.util.Set;

import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.project.RunOutput;

/**
 * Class that provides some Diffany-specific assertion methods for JUnit tests
 * 
 * @author Sofie Van Landeghem
 */
public class TestGeneric
{
	
	/**
	 * Protected method that asserts whether a certain edge is present in a network.
	 * This method will fail during JUnit testing when there is no edge or more than one edge between the specified node names.
	 * This method will also fail if the found edge has the wrong symmetry, weight, negation, or interaction type.
	 * 
	 * @param n the network to test
	 * @param sourceID the ID of the source node in the network
	 * @param targetID the ID of the target node in the network
	 * @param symm whether or not the relationship is symmetrical
	 * @param type the type the edge should have
	 * @param negated whether or not the edge should be negated
	 * @param weight the weight the edge should have
	 */
	protected void assertAnEdge(Network n, String sourceID, String targetID, boolean symm, String type, boolean negated, double weight)
	{
		Set<Edge> edges = n.getAllEdges(sourceID, targetID);
		boolean found = false;
		for (Edge edge : edges)
		{
			if (edge.isSymmetrical() == symm && edge.isNegated() == negated)
			{
				if (edge.getType().equals(type))
				{
					if ((edge.getWeight() < weight + 0.00000005) && (edge.getWeight() > weight - 0.00000005))
					{
						found = true;
					}
				}
			}
		}
		assertEquals(found, true);
	}
	
	/**
	 * Protected method that asserts the number of differential output pairs in the output result (may be 0).
	 * @param output the output of the run
	 * @param number the expected number of output pairs
	 */
	protected void assertNrPairs(RunOutput output, int number)
	{
		int pairs = output.getOutputAsPairs().size();
		assertEquals(number, pairs);
	}
	
	/**
	 * Protected method that asserts the number of differential networks in the output result (may be 0).
	 * @param output the output of the run
	 * @param number the expected number of differential networks
	 */
	protected void assertNrDiffNetworks(RunOutput output, int number)
	{
		int diffs = output.getDifferentialNetworks().size();
		assertEquals(number, diffs);
	}
	
	
	/**
	 * Protected method that asserts the number of consensus networks in the output result (may be 0).
	 * @param output the output of the run
	 * @param number the expected number of consensus networks
	 */
	protected void assertNrConsensusNetworks(RunOutput output, int number)
	{
		int consensusNr = output.getConsensusNetworks().size();
		assertEquals(number, consensusNr);
	}

}
