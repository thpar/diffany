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
import be.svlandeg.diffany.examples.FuzzyNetworks;

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
	public void testFuzzyOverlapWithReferenceAsCondition_4()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 4;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff, true);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 15, true);

		// Testing that there is exactly one overlap network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "X", "Y", false, "regulation", false, 0.5);
	}
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 3 out of 4 (75%).
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyOverlapWithReferenceAsCondition_3()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 3;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff, true);

		
		boolean exception = false;
		try
		{
			new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 10, true);		// this is an ID of an input network and should thus not work!
		}
		catch (IllegalArgumentException e)
		{
			exception = true;
		}
		assertTrue(exception);
		
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 20, true);	

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(3, sEdges.size());

		assertAnEdge(on, "A", "B", false, "colocalization", false, 0.3);
		assertAnEdge(on, "X", "Y", false, "regulation", false, 0.6);
		assertAnEdge(on, "M", "N", false, "phosphorylation", true, 0.3);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 3 out of 4 (75%).
	 * Instead of the minimum operator, the maximum is used.
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyOverlapWithReferenceAsCondition_3_maxOperator()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 3;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff, true);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 30, false);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(3, sEdges.size());

		assertAnEdge(on, "A", "B", false, "colocalization", false, 0.8);
		assertAnEdge(on, "X", "Y", false, "regulation", false, 0.8);
		assertAnEdge(on, "M", "N", false, "phosphorylation", true, 0.6);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 2 out of 4 (50%).
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyOverlapWithReferenceAsCondition_2()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 2;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff, true);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 40, true);
		
		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(4, sEdges.size());

		assertAnEdge(on, "A", "B", false, "colocalization", false, 0.6);
		assertAnEdge(on, "B", "A", false, "ppi", false, 0.4);

		assertAnEdge(on, "X", "Y", false, "positive_regulation", false, 0.7);	
		
		assertAnEdge(on, "M", "N", false, "ptm", true, 0.5);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 1 out of 4 (25%).
	 * This should throw an error, as the overlap factor should always be above 1 (i.e. at least 2).
	 */
	@Test
	public void testFuzzyOverlapWithReferenceAsCondition_1()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		Project p = ex.getProject();
		int overlap_cutoff = 1;
		
		boolean exception = false;
		try 
		{
			ex.getTestConfigurationWithoutReference(p, overlap_cutoff, true);
		}
		catch(IllegalArgumentException e)
		{
			exception = true;
		}
		assertTrue(exception);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 0 out of 4 (0%). 
	 * This should throw an error, as the overlap factor should always be above 1 (i.e. at least 2).
	 */
	@Test
	public void testFuzzyOverlapWithReferenceAsCondition_0()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		Project p = ex.getProject();
		int overlap_cutoff = 0;

		boolean exception = false;
		try
		{
			ex.getTestConfigurationWithoutReference(p, overlap_cutoff, true);
		}
		catch (IllegalArgumentException e)
		{
			exception = true;
		}
		assertTrue(exception);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 4.
	 * This should throw an error, as there are only 3 input networks.
	 */
	@Test
	public void testFuzzyOverlapWithoutReference_4()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 4;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff, false);

		boolean exception = false;
		try
		{
			new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 15, true);
		}
		catch (IllegalArgumentException e)
		{
			exception = true;
		}
		assertTrue(exception);
	}
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 3 out of 3 (100%).
	 * This method uses only the condition-specific networks, all cast as generic input networks.
	 */
	@Test
	public void testFuzzyOverlapWithoutReference_3()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 3;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff, false);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 20, true);	

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "X", "Y", false, "regulation", false, 0.6);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 3 out of 3 (100%).
	 * Instead of the minimum operator, the maximum is used.
	 * This method uses only the condition-specific networks, all cast as generic input networks.
	 */
	@Test
	public void testFuzzyOverlapWithoutReference_3_maxOperator()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 3;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff, false);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 30, false);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "X", "Y", false, "regulation", false, 0.8);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 2 out of 3 (67%).
	 * This method uses only the condition-specific networks, all cast as generic input networks.
	 */
	@Test
	public void testFuzzyOverlapWithoutReference_2()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 2;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff, false);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 40, true);
		
		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(4, sEdges.size());

		assertAnEdge(on, "A", "B", false, "ppi", false, 0.3);
		assertAnEdge(on, "B", "A", false, "ppi", false, 0.4);

		assertAnEdge(on, "X", "Y", false, "positive_regulation", false, 0.7);	
		
		assertAnEdge(on, "M", "N", false, "phosphorylation", true, 0.3);
	}


	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results when varying the fuzziness overlap factor.
	 * This method defines one of the networks to be a reference network.
	 */
	@Test
	public void testFuzzyOverlapWithReference_4()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 4;
		int ID = ex.getTestConfigurationWithReference(p, overlap_cutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 70, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "X", "Y", false, "regulation", false, 0.5);
	}
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results when varying the fuzziness overlap factor.
	 * This method defines one of the networks to be a reference network.
	 */
	@Test
	public void testFuzzyOverlapWithReference_3()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 3;
		int ID = ex.getTestConfigurationWithReference(p, overlap_cutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 70, true);

		// Testing that there is exactly one overlap network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(3, sEdges.size());

		assertAnEdge(on, "A", "B", false, "colocalization", false, 0.3);
		assertAnEdge(on, "X", "Y", false, "regulation", false, 0.5);		// important different with scenario without reference (0.6) !!!
		assertAnEdge(on, "M", "N", false, "phosphorylation", true, 0.3);
	}

	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results when varying the fuzziness overlap factor.
	 * This method defines one of the networks to be a reference network.
	 */
	@Test
	public void testFuzzyOverlapWithReference_2()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 2;
		int ID = ex.getTestConfigurationWithReference(p, overlap_cutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 70, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(3, sEdges.size());

		assertAnEdge(on, "A", "B", false, "colocalization", false, 0.6);
		assertAnEdge(on, "X", "Y", false, "negative_regulation", false, 0.5);		// important different with scenario without reference (pos 0.7) !!!
		assertAnEdge(on, "M", "N", false, "ptm", true, 0.5);
	}

	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results when varying the fuzziness overlap factor.
	 * This should throw an error, as the overlap factor should always be above 1 (i.e. at least 2).
	 */
	@Test
	public void testFuzzyOverlapWithReference_1()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		Project p = ex.getProject();
		int overlap_cutoff = 1;
		boolean exception = false;
		try
		{
			ex.getTestConfigurationWithReference(p, overlap_cutoff);
		}
		catch (IllegalArgumentException e)
		{
			exception = true;
		}
		assertTrue(exception);
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

	/**
	 * Private method that asserts the number of overlap networks in the output result (may be 0).
	 */
	private void assertNrOverlapNetworks(RunOutput output, int number)
	{
		int overlaps = output.getOverlappingNetworks().size();
		assertEquals(number, overlaps);
	}
}
