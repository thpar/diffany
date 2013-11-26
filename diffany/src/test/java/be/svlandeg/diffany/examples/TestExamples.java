package be.svlandeg.diffany.examples;

import static org.junit.Assert.assertEquals;

import java.util.Collection;
import java.util.Set;

import org.junit.Test;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.*;

/** 
 * Class that automatically tests the outputs of the small examples 
 * contained in the be.svlandeg.diffany.examples package.
 * 
 * @author Sofie Van Landeghem
 */
public class TestExamples 
{
	
	
	/**
	 * Test whether the small example figure 1C from the Bandyopadhyay et al. 2010 paper
	 * produces correct results.
	 */
	@Test
	public void testBandyopadhyay()
	{
		Bandyopadhyay2010 ex = new Bandyopadhyay2010();
		double cutoff = 0.25;
		Project p = ex.getProjectFigure1C();
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, cutoff);
		
		// Testing that there is exactly one differential network created
		Collection<DifferentialNetwork> dNetworks = p.getDifferentialNetworks();
		assertEquals(1, dNetworks.size());
		
		// Testing the edges in the differential network
		DifferentialNetwork dNetwork = dNetworks.iterator().next();
		Set<Edge> dEdges =  dNetwork.getEdges();
		assertEquals(3, dEdges.size());
		
		assertOneEdge(dNetwork, "A", "B", true, false, "increase", false, 0.7);
		assertOneEdge(dNetwork, "A", "C", true, false, "decrease", false, 0.7);
		assertOneEdge(dNetwork, "C", "E", true, false, "decrease", false, 0.8);
		
		// Testing the edges in the corresponding shared network
		SharedNetwork sNetwork = dNetwork.getSharedNetwork();
		Set<Edge> sEdges =  sNetwork.getEdges();
		assertEquals(1, sEdges.size());
		
		assertOneEdge(sNetwork, "A", "D", true, false, "negative", false, 0.9);
	}
	
	/**
	 * Test whether the small example figure 3A from the Ideker et al. 2011 paper
	 * produces correct results.
	 */
	@Test
	public void testIdeker()
	{
		Ideker2011 ex = new Ideker2011();
		double cutoff = 0.25;
		Project p = ex.getProjectFigure3A();
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, cutoff);
		
		// Testing that there is exactly one differential network created
		Collection<DifferentialNetwork> dNetworks = p.getDifferentialNetworks();
		assertEquals(1, dNetworks.size());
		
		// Testing the edges in the differential network
		DifferentialNetwork dNetwork = dNetworks.iterator().next();
		Set<Edge> dEdges =  dNetwork.getEdges();
		assertEquals(3, dEdges.size());
		
		assertOneEdge(dNetwork, "A", "B", true, false, "increase", false, 0.7);
		assertOneEdge(dNetwork, "A", "C", true, false, "decrease", false, 1.2);
		assertOneEdge(dNetwork, "A", "E", true, false, "decrease", false, 0.8);
		
		// Testing the edges in the corresponding shared network
		SharedNetwork sNetwork = dNetwork.getSharedNetwork();
		Set<Edge> sEdges =  sNetwork.getEdges();
		assertEquals(2, sEdges.size());
		
		assertOneEdge(sNetwork, "A", "D", true, false, "negative", false, 0.7);
		assertOneEdge(sNetwork, "A", "F", true, false, "negative", false, 1);
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
	private void assertOneEdge(Network n, String sourceName, String targetName, boolean symm, boolean normalized, String type, boolean negated, double weight)
	{
		Set<Edge> edges = n.getAllEdgesByName(sourceName, targetName, symm, normalized);
		assertEquals(1, edges.size());
		Edge edge = edges.iterator().next();
		assertEquals(edge.getType(), type);
		assertEquals(edge.isSymmetrical(), symm);
		assertEquals(edge.isNegated(), negated);
		assertEquals(edge.getWeight(), weight, 0);
	}
}
