package be.svlandeg.diffany.junit;

import static org.junit.Assert.*;

import java.util.Set;

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.OverlappingNetwork;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.examples.FuzzyOverlap;

/**
 * Class that automatically tests the outputs of some examples of 'fuzzy' overlap networks.
 * 
 * @author Sofie Van Landeghem
 */
public class TestFuzzyOverlap
{

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 4 out of 4 100%.
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyOverlapWithoutReference_4()
	{
		FuzzyOverlap ex = new FuzzyOverlap();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 4;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, false, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		/*for (Edge e : sEdges)
		{
			System.out.println(e);
		}*/
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "X", "Y", false, false, "regulation", false, 0.3);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 3 out of 4 (75%).
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyOverlapWithoutReference_3()
	{
		FuzzyOverlap ex = new FuzzyOverlap();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 3;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, false, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		/*for (Edge e : sEdges)
		{
			System.out.println(e);
		}*/
		assertEquals(2, sEdges.size());

		assertAnEdge(on, "A", "B", false, false, "ppi", false, 0.3);
		assertAnEdge(on, "X", "Y", false, false, "regulation", false, 0.3);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 2 out of 4 (50%).
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyOverlapWithoutReference_2()
	{
		FuzzyOverlap ex = new FuzzyOverlap();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 2;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, false, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		/*for (Edge e : sEdges)
		{
			System.out.println(e);
		}*/
		assertEquals(3, sEdges.size());

		assertAnEdge(on, "A", "B", false, false, "ppi", false, 0.6);
		assertAnEdge(on, "B", "A", false, false, "ppi", false, 0.4);

		//assertAnEdge(on, "X", "Y", false, false, "positive regulation", false, 0.3);	will be removed after clean-up of the network because 'regulation' has higher weight
		assertAnEdge(on, "X", "Y", false, false, "regulation", false, 0.6);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 1 out of 4 (25%).
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyOverlapWithoutReference_1()
	{
		FuzzyOverlap ex = new FuzzyOverlap();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 1;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, false, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		/*for (Edge e : sEdges)
		{
			System.out.println(e);
		}*/
		assertEquals(3, sEdges.size());

		assertAnEdge(on, "A", "B", false, false, "ppi", false, 0.8);
		assertAnEdge(on, "B", "A", false, false, "ppi", false, 0.8);

		assertAnEdge(on, "X", "Y", false, false, "positive regulation", false, 0.8);

		// current cleaning will take only the edge with the highest weight!
		//assertAnEdge(on, "X", "Y", false, false, "negative regulation", false, 0.5);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 0 out of 4 (0%). This should not give any results.
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyOverlapWithoutReference_0()
	{
		FuzzyOverlap ex = new FuzzyOverlap();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 0;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff);
		boolean exception = false;

		try
		{
			new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, false, true);
		}
		catch (IllegalArgumentException e)
		{
			exception = true;
		}
		assertTrue(exception);

		// Testing that there is no differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 0);
	}

	/**
	 * TODO: different factors
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results when varying the fuzziness overlap factor.
	 * This method defines one of the networks to be a reference network.
	 */
	@Test
	public void testFuzzyOverlapWithReference()
	{
		FuzzyOverlap ex = new FuzzyOverlap();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 4;
		int ID = ex.getTestConfigurationWithReference(p, overlap_cutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, false, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "X", "Y", false, false, "regulation", false, 0.3);
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
	private void assertAnEdge(Network n, String sourceName, String targetName, boolean symm, boolean normalized, String type, boolean negated, double weight)
	{
		Set<Edge> edges = n.getAllEdgesByName(sourceName.toLowerCase(), targetName.toLowerCase(), symm, normalized);
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
	 * Private method that asserts the number of overlap networks in the output result (may be 0).
	 */
	private void assertNrOverlapNetworks(RunOutput output, int number)
	{
		int overlaps = output.getOverlappingNetworks().size();
		assertEquals(number, overlaps);
	}
}
