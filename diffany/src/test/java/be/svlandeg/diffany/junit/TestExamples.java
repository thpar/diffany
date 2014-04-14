package be.svlandeg.diffany.junit;

import static org.junit.Assert.assertEquals;

import java.util.Collection;
import java.util.Set;

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.*;
import be.svlandeg.diffany.core.project.DifferentialOutput;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.examples.ActivityFlowTest;
import be.svlandeg.diffany.examples.Bandyopadhyay2010;
import be.svlandeg.diffany.examples.ConflictingEdgesTest;
import be.svlandeg.diffany.examples.Ideker2011;
import be.svlandeg.diffany.examples.MultipleConditionTest;
import be.svlandeg.diffany.examples.ProcessTest;

/**
 * Class that automatically tests the outputs of the small examples
 * contained in the be.svlandeg.diffany.junit package.
 * 
 * @author Sofie Van Landeghem
 */
public class TestExamples
{

	/**
	 * JUNIT Test: check whether the small example figure 1C from the Bandyopadhyay et al. 2010 paper
	 * produces correct results.
	 */
	@Test
	public void testBandyopadhyay()
	{
		Bandyopadhyay2010 ex = new Bandyopadhyay2010();
		double cutoff = 0.0;
		Project p = ex.getProjectFigure1C();
		int ID = ex.getTestConfiguration1C(p);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, cutoff, true, true);

		// Testing that there is exactly one differential network created
		Collection<DifferentialOutput> dNetworks = p.getRunConfiguration(ID).getDifferentialOutputs();
		assertEquals(1, dNetworks.size());

		// Testing the edges in the differential network
		DifferentialOutput output = dNetworks.iterator().next();
		OutputNetworkPair pair = output.getOutputAsPair();
		DifferentialNetwork dNetwork = pair.getDifferentialNetwork();
		Set<Edge> dEdges = dNetwork.getEdges();
		assertEquals(3, dEdges.size());

		assertAnEdge(dNetwork, "A", "B", true, false, "increase_genetic_interaction", false, 0.7);
		assertAnEdge(dNetwork, "A", "C", true, false, "decrease_genetic_interaction", false, 0.7);
		assertAnEdge(dNetwork, "C", "E", true, false, "decrease_genetic_interaction", false, 0.8);

		// Testing the edges in the corresponding overlapping network
		OverlappingNetwork sNetwork = pair.getOverlappingNetwork();
		Set<Edge> sEdges = sNetwork.getEdges();
		assertEquals(2, sEdges.size());

		assertAnEdge(sNetwork, "A", "D", true, false, "negative_genetic_interaction", false, 1.1);
		assertAnEdge(sNetwork, "A", "B", true, false, "genetic_interaction", false, 0.3);
	}

	/**
	 * JUNIT Test: check whether the small example figure 3A from the Ideker et al. 2011 paper
	 * produces correct results.
	 */
	@Test
	public void testIdeker()
	{
		Ideker2011 ex = new Ideker2011();
		double cutoff = 0.0;
		Project p = ex.getProjectFigure3A();
		int ID = ex.getTestConfiguration3A(p);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, cutoff, true, true);

		// Testing that there is exactly one differential network created
		Collection<DifferentialOutput> dNetworks = p.getRunConfiguration(ID).getDifferentialOutputs();
		assertEquals(1, dNetworks.size());

		// Testing the edges in the differential network
		DifferentialOutput output = dNetworks.iterator().next();
		OutputNetworkPair pair = output.getOutputAsPair();
		DifferentialNetwork dNetwork = pair.getDifferentialNetwork();

		Set<Edge> dEdges = dNetwork.getEdges();

		assertEquals(3, dEdges.size());

		assertAnEdge(dNetwork, "A", "B", true, false, "increase_genetic_interaction", false, 0.7);
		assertAnEdge(dNetwork, "A", "C", true, false, "decrease_genetic_interaction", false, 1.2);
		assertAnEdge(dNetwork, "A", "E", true, false, "decrease_genetic_interaction", false, 0.8);

		// Testing the edges in the corresponding overlapping network
		OverlappingNetwork sNetwork = pair.getOverlappingNetwork();
		Set<Edge> sEdges = sNetwork.getEdges();
		assertEquals(3, sEdges.size());

		assertAnEdge(sNetwork, "A", "D", true, false, "negative_genetic_interaction", false, 0.7);
		assertAnEdge(sNetwork, "A", "F", true, false, "negative_genetic_interaction", false, 1);
		assertAnEdge(sNetwork, "A", "B", true, false, "genetic_interaction", false, 0.3);
	}

	/**
	 * JUNIT Test: check whether the example activity flow network produces correct results.
	 */
	@Test
	public void testActivityFlowNetwork()
	{
		ActivityFlowTest ex = new ActivityFlowTest();
		double cutoff = 0.0;
		Project p = ex.getTestProject();
		int ID = ex.getTestConfiguration(p);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, cutoff, true, true);

		// Testing that there is exactly one differential network created
		Collection<DifferentialOutput> dNetworks = p.getRunConfiguration(ID).getDifferentialOutputs();
		assertEquals(1, dNetworks.size());

		// Testing the edges in the differential network
		DifferentialOutput output = dNetworks.iterator().next();
		OutputNetworkPair pair = output.getOutputAsPair();
		DifferentialNetwork dNetwork = pair.getDifferentialNetwork();

		Set<Edge> dEdges = dNetwork.getEdges();
		assertEquals(7, dEdges.size());

		assertAnEdge(dNetwork, "S", "T", false, false, "increases_regulation", false, 1);
		assertAnEdge(dNetwork, "K", "J", false, false, "increases_regulation", false, 2);
		assertAnEdge(dNetwork, "J", "K", false, false, "increases_regulation", false, 2);
		assertAnEdge(dNetwork, "A", "B", false, false, "increases_regulation", false, 1);
		assertAnEdge(dNetwork, "B", "A", false, false, "decreases_regulation", false, 1);
		assertAnEdge(dNetwork, "M", "N", false, false, "increases_regulation", false, 12);
		assertAnEdge(dNetwork, "N", "M", false, false, "increases_regulation", false, 7);

		// Testing the edges in the corresponding overlapping network
		OverlappingNetwork sNetwork = pair.getOverlappingNetwork();
		Set<Edge> sEdges = sNetwork.getEdges();
		assertEquals(5, sEdges.size());

		assertAnEdge(sNetwork, "A", "B", false, false, "positive_regulation", false, 2);
		assertAnEdge(sNetwork, "G", "H", false, false, "negative_regulation", true, 3);
		assertAnEdge(sNetwork, "H", "G", false, false, "negative_regulation", true, 3);
		assertAnEdge(sNetwork, "X", "Y", false, false, "positive_regulation", true, 2);
		assertAnEdge(sNetwork, "M", "N", false, false, "regulation", false, 5);
	}

	/**
	 * JUNIT Test: check whether the example process network produces correct results.
	 */
	@Test
	public void testProcessNetwork()
	{
		ProcessTest ex = new ProcessTest();
		double cutoff = 0.0;
		Project p = ex.getTestProject();
		int ID = ex.getTestConfiguration(p);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, cutoff, true, true);

		// Testing that there is exactly one differential network created
		Collection<DifferentialOutput> dNetworks = p.getRunConfiguration(ID).getDifferentialOutputs();
		assertEquals(1, dNetworks.size());

		// Testing the edges in the differential network
		DifferentialOutput output = dNetworks.iterator().next();
		OutputNetworkPair pair = output.getOutputAsPair();
		DifferentialNetwork dNetwork = pair.getDifferentialNetwork();

		Set<Edge> dEdges = dNetwork.getEdges();

		assertEquals(9, dEdges.size());

		assertAnEdge(dNetwork, "X", "Y", false, false, "decreases_ptm", false, 1);
		assertAnEdge(dNetwork, "Y", "X", false, false, "decreases_ptm", false, 1);
		assertAnEdge(dNetwork, "A", "B", true, false, "decrease_ppi", false, 2);
		assertAnEdge(dNetwork, "G", "H", false, false, "increases_ptm", false, 4);
		assertAnEdge(dNetwork, "H", "G", false, false, "decreases_ubiquitination", false, 1);
		assertAnEdge(dNetwork, "M", "N", true, false, "increase_ppi", false, 3);
		assertAnEdge(dNetwork, "S", "T", true, false, "decrease_ppi", false, 3);
		assertAnEdge(dNetwork, "S", "T", false, false, "increases_phosphorylation", false, 2);
		assertAnEdge(dNetwork, "T", "S", false, false, "increases_phosphorylation", false, 2);

		// Testing the edges in the corresponding overlapping network
		OverlappingNetwork sNetwork = pair.getOverlappingNetwork();
		Set<Edge> sEdges = sNetwork.getEdges();
		assertEquals(5, sEdges.size());

		assertAnEdge(sNetwork, "X", "Y", false, false, "ptm", false, 3);
		assertAnEdge(sNetwork, "Y", "X", false, false, "ptm", false, 3);
		assertAnEdge(sNetwork, "G", "H", false, false, "ptm", false, 1);
		assertAnEdge(sNetwork, "K", "J", false, false, "phosphorylation", true, 4);
		assertAnEdge(sNetwork, "J", "K", false, false, "phosphorylation", true, 4);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results.
	 */
	@Test
	public void testMultipleConditions()
	{
		MultipleConditionTest ex = new MultipleConditionTest();
		double cutoff = 0.0;
		Project p = ex.getTestProject();
		int ID = ex.getTestDiffConfiguration(p);
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, cutoff, true, true);

		// Testing that there is exactly one differential network created
		Collection<DifferentialOutput> dNetworks = p.getRunConfiguration(ID).getDifferentialOutputs();
		assertEquals(1, dNetworks.size());

		// Testing the edges in the differential network
		DifferentialOutput output = dNetworks.iterator().next();
		OutputNetworkPair pair = output.getOutputAsPair();
		DifferentialNetwork dNetwork = pair.getDifferentialNetwork();

		Set<Edge> dEdges = dNetwork.getEdges();

		assertEquals(8, dEdges.size());

		assertAnEdge(dNetwork, "W", "Z", true, false, "decrease_ppi", false, 0.5);
		assertAnEdge(dNetwork, "A", "B", true, false, "decrease_ppi", false, 0.3);
		assertAnEdge(dNetwork, "A", "D", true, false, "increase_ppi", false, 0.75);
		assertAnEdge(dNetwork, "A", "Z", true, false, "decrease_ppi", false, 0.8);

		assertAnEdge(dNetwork, "A", "B", false, false, "decreases_phosphorylation", false, 1.1);

		assertAnEdge(dNetwork, "M", "N", false, false, "increases_phosphorylation", false, 4);
		assertAnEdge(dNetwork, "P", "M", false, false, "increases_ptm", false, 2);
		assertAnEdge(dNetwork, "N", "P", false, false, "decreases_phosphorylation", false, 3);

		// Testing the edges in the corresponding overlapping network
		OverlappingNetwork sNetwork = pair.getOverlappingNetwork();
		Set<Edge> sEdges = sNetwork.getEdges();
		assertEquals(3, sEdges.size());

		assertAnEdge(sNetwork, "A", "B", true, false, "ppi", false, 0.3);
		assertAnEdge(sNetwork, "A", "C", true, false, "ppi", false, 0.6);
		assertAnEdge(sNetwork, "M", "N", false, false, "phosphorylation", false, 2);
	}

	/**
	 * JUNIT Test: check whether the example network with edge conflicts produces correct results.
	 */
	@Test
	public void testConflicts()
	{
		ConflictingEdgesTest ex = new ConflictingEdgesTest();
		double cutoff = 0.0;
		Project p = ex.getTestProject();
		int ID = ex.getTestConfiguration(p);
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, cutoff, true, true);

		// Testing that there is exactly one differential network created
		Collection<DifferentialOutput> dNetworks = p.getRunConfiguration(ID).getDifferentialOutputs();
		assertEquals(1, dNetworks.size());

		// Testing the edges in the differential network
		DifferentialOutput output = dNetworks.iterator().next();
		OutputNetworkPair pair = output.getOutputAsPair();
		DifferentialNetwork dNetwork = pair.getDifferentialNetwork();

		Set<Edge> dEdges = dNetwork.getEdges();

		assertEquals(9, dEdges.size());

		assertAnEdge(dNetwork, "A", "B", false, false, "increases_regulation", false, 6);
		assertAnEdge(dNetwork, "A", "B", false, false, "decreases_ptm", false, 5);
		assertAnEdge(dNetwork, "A", "B", false, false, "decreases_somerandomInteraction", false, 4);
		assertAnEdge(dNetwork, "G", "H", false, false, "decreases_unspecified_regulation", false, 3);
		assertAnEdge(dNetwork, "G", "H", false, false, "decreases_ptm", false, 1);
		assertAnEdge(dNetwork, "J", "K", false, false, "decreases_unspecified_regulation", false, 3);
		assertAnEdge(dNetwork, "K", "J", false, false, "decreases_regulation", false, 4);
		assertAnEdge(dNetwork, "J", "K", false, false, "increases_ptm", false, 6);
		assertAnEdge(dNetwork, "K", "J", false, false, "increases_ptm", false, 1);

		// Testing the edges in the corresponding overlapping network
		OverlappingNetwork sNetwork = pair.getOverlappingNetwork();
		Set<Edge> sEdges = sNetwork.getEdges();
		assertEquals(5, sEdges.size());

		assertAnEdge(sNetwork, "A", "B", false, false, "positive_regulation", false, 2);
		assertAnEdge(sNetwork, "G", "H", false, false, "regulation", false, 4);
		assertAnEdge(sNetwork, "G", "H", false, false, "ptm", false, 2);
		assertAnEdge(sNetwork, "J", "K", false, false, "regulation", false, 4);
		assertAnEdge(sNetwork, "K", "J", false, false, "ptm", false, 2);
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
		Set<Edge> edges = n.getAllEdgesByName(sourceName, targetName, symm, normalized);
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
