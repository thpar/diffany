package be.svlandeg.diffany.junit;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.examples.FuzzyNetworks;


/**
 * Class that automatically tests the outputs of some examples of 'fuzzy' differential networks.
 * 
 * @author Sofie Van Landeghem
 */
public class TestFuzzyDiff
{
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results when varying the fuzziness overlap factor.
	 * This method defines one of the networks to be a reference network.
	 */
	@Test
	public void testFuzzyDiff_3()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 3;
		int ID = ex.getTestConfigurationWithReference(p, overlap_cutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, 70, -1, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 0);
		assertNrDifferentialNetworks(output, 1);

		// Testing the edges in the overlap network
		DifferentialNetwork dn = output.getDifferentialNetworks().iterator().next();

		Set<Edge> sEdges = dn.getEdges();
		assertEquals(1, sEdges.size());

		assertAnEdge(dn, "X", "Y", false, "increases_unspecified_regulation", false, 0.1);
	}
	
	/**
	 * Private method that asserts the number of overlap networks in the output result (may be 0).
	 */
	private void assertNrOverlapNetworks(RunOutput output, int number)
	{
		int overlaps = output.getOverlappingNetworks().size();
		assertEquals(number, overlaps);
	}
	
	/**
	 * Private method that asserts the number of overlap networks in the output result (may be 0).
	 */
	private void assertNrDifferentialNetworks(RunOutput output, int number)
	{
		int diffs = output.getDifferentialNetworks().size();
		assertEquals(number, diffs);
	}
	
	/**
	 * Private method that asserts whether a certain edge is present in a network.
	 * This method will fail during JUnit testing when there is no edge or more than one edge between the specified node names.
	 * This method will also fail if the found edge has the wrong symmetry, weight, negation, or interaction type.
	 * 
	 * @param n the network to test
	 * @param sourceName the name of the source node in the network
	 * @param targetName the name of the target node in the network
	 * @param symm whether or not the relationship is symmetrical
	 * @param normalized whether or not the node names are in normalized form
	 * @param type the type the edge should have
	 * @param negated whether or not the edge should be negated
	 * @param weight the weight the edge should have
	 */
	private void assertAnEdge(Network n, String sourceName, String targetName, boolean symm, String type, boolean negated, double weight)
	{
		Set<Edge> edges = n.getAllEdgesByName(sourceName.toLowerCase(), targetName.toLowerCase(), symm);
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

}
