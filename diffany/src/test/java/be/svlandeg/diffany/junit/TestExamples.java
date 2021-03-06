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

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.progress.ProgressListener;
import be.svlandeg.diffany.core.progress.StandardProgressListener;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunOutput;
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
public class TestExamples extends TestGeneric
{

	/**
	 * JUNIT Test: check whether the small example figure 1C from the Bandyopadhyay et al. 2010 paper
	 * produces correct results.
	 */
	@Test
	public void testBandyopadhyay()
	{
		Bandyopadhyay2010 ex = new Bandyopadhyay2010();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int ID = ex.getDefaultRunConfigurationID(p);
		ProgressListener listener = new StandardProgressListener(false);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, weight_cutoff, true, true, 10, true, listener);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 1);

		// Testing the edges in the differential network
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
		DifferentialNetwork dNetwork = pair.getDifferentialNetwork();
		Set<Edge> dEdges = dNetwork.getEdges();
		assertEquals(3, dEdges.size());

		assertAnEdge(dNetwork, "A", "B", true, "increase_genetic_interaction", false, 0.7);
		assertAnEdge(dNetwork, "A", "C", true, "decrease_genetic_interaction", false, 0.7);
		assertAnEdge(dNetwork, "C", "E", true, "decrease_genetic_interaction", false, 0.8);

		// Testing the edges in the corresponding consensus network
		ConsensusNetwork sNetwork = pair.getConsensusNetwork();
		Set<Edge> sEdges = sNetwork.getEdges();
		assertEquals(2, sEdges.size());

		assertAnEdge(sNetwork, "A", "D", true, "negative_genetic_interaction", false, 1.1);
		assertAnEdge(sNetwork, "A", "B", true, "genetic_interaction", false, 0.3);
	}

	/**
	 * JUNIT Test: check whether the small example figure 3A from the Ideker et al. 2011 paper
	 * produces correct results.
	 */
	@Test
	public void testIdeker()
	{
		Ideker2011 ex = new Ideker2011();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int ID = ex.getDefaultRunConfigurationID(p);
		ProgressListener listener = new StandardProgressListener(false);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, weight_cutoff, true, true, 10, true, listener);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 1);

		// Testing the edges in the differential network
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
		DifferentialNetwork dNetwork = pair.getDifferentialNetwork();

		Set<Edge> dEdges = dNetwork.getEdges();

		assertEquals(3, dEdges.size());

		assertAnEdge(dNetwork, "A", "B", true, "increase_genetic_interaction", false, 0.7);
		assertAnEdge(dNetwork, "A", "C", true, "decrease_genetic_interaction", false, 1.2);
		assertAnEdge(dNetwork, "A", "E", true, "decrease_genetic_interaction", false, 0.8);

		// Testing the edges in the corresponding consensus network
		ConsensusNetwork sNetwork = pair.getConsensusNetwork();
		Set<Edge> sEdges = sNetwork.getEdges();
		assertEquals(3, sEdges.size());

		assertAnEdge(sNetwork, "A", "D", true, "negative_genetic_interaction", false, 0.7);
		assertAnEdge(sNetwork, "A", "F", true, "negative_genetic_interaction", false, 1);
		assertAnEdge(sNetwork, "A", "B", true, "genetic_interaction", false, 0.3);
	}

	/**
	 * JUNIT Test: check whether the example activity flow network produces correct results.
	 */
	@Test
	public void testActivityFlowNetwork()
	{
		ActivityFlowTest ex = new ActivityFlowTest();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int ID = ex.getDefaultRunConfigurationID(p);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, weight_cutoff, true, true, 10, true, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 1);

		// Testing the edges in the differential network
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
		DifferentialNetwork dNetwork = pair.getDifferentialNetwork();

		Set<Edge> dEdges = dNetwork.getEdges();

		assertEquals(7, dEdges.size());

		assertAnEdge(dNetwork, "S", "T", false, "increases_regulation", false, 1);
		assertAnEdge(dNetwork, "K", "J", false, "increases_regulation", false, 2);
		assertAnEdge(dNetwork, "J", "K", false, "increases_regulation", false, 2);
		assertAnEdge(dNetwork, "A", "B", false, "increases_regulation", false, 1);
		assertAnEdge(dNetwork, "B", "A", false, "decreases_regulation", false, 1);
		assertAnEdge(dNetwork, "M", "N", false, "increases_regulation", false, 12);
		assertAnEdge(dNetwork, "N", "M", false, "increases_regulation", false, 7);

		// Testing the edges in the corresponding consensus network
		ConsensusNetwork sNetwork = pair.getConsensusNetwork();
		Set<Edge> sEdges = sNetwork.getEdges();
		assertEquals(5, sEdges.size());

		assertAnEdge(sNetwork, "A", "B", false, "positive_regulation", false, 2);
		assertAnEdge(sNetwork, "G", "H", false, "negative_regulation", true, 3);
		assertAnEdge(sNetwork, "H", "G", false, "negative_regulation", true, 3);
		assertAnEdge(sNetwork, "X", "Y", false, "positive_regulation", true, 2);
		assertAnEdge(sNetwork, "M", "N", false, "regulation", false, 5);
	}

	/**
	 * JUNIT Test: check whether the example process network produces correct results.
	 */
	@Test
	public void testProcessNetwork()
	{
		ProcessTest ex = new ProcessTest();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int ID = ex.getDefaultRunConfigurationID(p);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, weight_cutoff, true, true, 10, true, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 1);

		// Testing the edges in the differential network
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
		DifferentialNetwork dNetwork = pair.getDifferentialNetwork();

		Set<Edge> dEdges = dNetwork.getEdges();

		assertEquals(9, dEdges.size());

		assertAnEdge(dNetwork, "X", "Y", false, "decreases_ptm", false, 1);
		assertAnEdge(dNetwork, "Y", "X", false, "decreases_ptm", false, 1);
		assertAnEdge(dNetwork, "A", "B", true, "decrease_ppi", false, 2);
		assertAnEdge(dNetwork, "G", "H", false, "increases_ptm", false, 4);
		assertAnEdge(dNetwork, "H", "G", false, "decreases_ubiquitination", false, 1);
		assertAnEdge(dNetwork, "M", "N", true, "increase_ppi", false, 3);
		assertAnEdge(dNetwork, "S", "T", true, "decrease_ppi", false, 3);
		assertAnEdge(dNetwork, "S", "T", false, "increases_phosphorylation", false, 2);
		assertAnEdge(dNetwork, "T", "S", false, "increases_phosphorylation", false, 2);

		// Testing the edges in the corresponding consensus network
		ConsensusNetwork sNetwork = pair.getConsensusNetwork();
		Set<Edge> sEdges = sNetwork.getEdges();
		assertEquals(5, sEdges.size());

		assertAnEdge(sNetwork, "X", "Y", false, "ptm", false, 3);
		assertAnEdge(sNetwork, "Y", "X", false, "ptm", false, 3);
		assertAnEdge(sNetwork, "G", "H", false, "ptm", false, 1);
		assertAnEdge(sNetwork, "K", "J", false, "phosphorylation", true, 4);
		assertAnEdge(sNetwork, "J", "K", false, "phosphorylation", true, 4);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results.
	 * This method specifically checks the 1-to-many algorithm.
	 */
	@Test
	public void testMultipleConditions1toMany()
	{
		MultipleConditionTest ex = new MultipleConditionTest();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int ID = ex.getDefaultRunConfigurationID(p);
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, null, null, 10, 20, true, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 1);

		// Testing the edges in the differential network
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
		DifferentialNetwork dNetwork = pair.getDifferentialNetwork();

		Set<Edge> dEdges = dNetwork.getEdges();

		assertEquals(8, dEdges.size());

		assertAnEdge(dNetwork, "W", "Z", true, "decrease_ppi", false, 0.5);
		assertAnEdge(dNetwork, "A", "B", true, "decrease_ppi", false, 0.3);
		assertAnEdge(dNetwork, "A", "D", true, "increase_ppi", false, 0.75);
		assertAnEdge(dNetwork, "A", "Z", true, "decrease_ppi", false, 0.8);

		assertAnEdge(dNetwork, "A", "B", false, "decreases_phosphorylation", false, 1.1);

		assertAnEdge(dNetwork, "M", "N", false, "increases_phosphorylation", false, 4);
		assertAnEdge(dNetwork, "P", "M", false, "increases_ptm", false, 2);
		assertAnEdge(dNetwork, "N", "P", false, "decreases_phosphorylation", false, 3);

		// Testing the edges in the corresponding consensus network
		ConsensusNetwork cNetwork = pair.getConsensusNetwork();
		Set<Edge> cEdges = cNetwork.getEdges();
		assertEquals(3, cEdges.size());

		assertAnEdge(cNetwork, "A", "B", true, "ppi", false, 0.3);
		assertAnEdge(cNetwork, "A", "C", true, "ppi", false, 0.6);
		assertAnEdge(cNetwork, "M", "N", false, "phosphorylation", false, 2);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results.
	 * This method specifically checks the pairwise algorithms of reference vs. conditions.
	 */
	@Test
	public void testMultipleConditionsPairwiseDiff()
	{
		MultipleConditionTest ex = new MultipleConditionTest();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int ID = ex.getDefaultRunConfigurationID(p);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, weight_cutoff, true, true, 10, true, null);

		// Testing that there are exactly two differential networks created (1 for each condition)

		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 2);

		// Testing the edges in the differential networks
		Map<String, OutputNetworkPair> outputs = new HashMap<String, OutputNetworkPair>();
		for (OutputNetworkPair pair : output.getOutputAsPairs())
		{
			outputs.put(pair.getDifferentialNetwork().getName(), pair);
		}

		// Salt vs. reference
		OutputNetworkPair saltPair = outputs.get("diff_Salty");
		DifferentialNetwork saltDiff = saltPair.getDifferentialNetwork();

		Set<Edge> dEdgesS = saltDiff.getEdges();

		assertEquals(11, dEdgesS.size());

		assertAnEdge(saltDiff, "W", "Z", true, "decrease_ppi", false, 0.5);
		assertAnEdge(saltDiff, "A", "Z", true, "decrease_ppi", false, 0.8);
		assertAnEdge(saltDiff, "A", "C", true, "decrease_ppi", false, 0.2);
		assertAnEdge(saltDiff, "A", "B", true, "decrease_ppi", false, 0.3);
		assertAnEdge(saltDiff, "A", "D", true, "increase_ppi", false, 0.9);
		assertAnEdge(saltDiff, "F", "D", true, "increase_ppi", false, 0.3);

		assertAnEdge(saltDiff, "A", "B", false, "decreases_phosphorylation", false, 1.1);

		assertAnEdge(saltDiff, "M", "N", false, "increases_phosphorylation", false, 6);
		assertAnEdge(saltDiff, "P", "M", false, "increases_phosphorylation", false, 2);
		assertAnEdge(saltDiff, "N", "P", false, "decreases_phosphorylation", false, 3);
		assertAnEdge(saltDiff, "O", "P", false, "increases_phosphorylation", false, 4);

		// Testing the edges in the corresponding consensus network
		ConsensusNetwork saltConsensus = saltPair.getConsensusNetwork();
		Set<Edge> cEdgesS = saltConsensus.getEdges();
		assertEquals(5, cEdgesS.size());

		assertAnEdge(saltConsensus, "A", "Z", true, "ppi", false, 0.1);
		assertAnEdge(saltConsensus, "A", "B", true, "ppi", false, 0.4);
		assertAnEdge(saltConsensus, "A", "C", true, "ppi", false, 0.6);

		assertAnEdge(saltConsensus, "M", "N", false, "phosphorylation", false, 2);
		assertAnEdge(saltConsensus, "M", "O", false, "phosphorylation", true, 1);

		// Draught vs. reference
		OutputNetworkPair draughtPair = outputs.get("diff_Draughty");
		DifferentialNetwork draughtDiff = draughtPair.getDifferentialNetwork();

		Set<Edge> dEdges = draughtDiff.getEdges();

		assertEquals(12, dEdges.size());

		assertAnEdge(draughtDiff, "W", "Z", true, "decrease_ppi", false, 0.5);
		assertAnEdge(draughtDiff, "A", "Z", true, "decrease_ppi", false, 0.9);
		assertAnEdge(draughtDiff, "A", "C", true, "increase_ppi", false, 0.4);
		assertAnEdge(draughtDiff, "A", "B", true, "decrease_ppi", false, 0.4);
		assertAnEdge(draughtDiff, "A", "D", true, "increase_ppi", false, 0.75);
		assertAnEdge(draughtDiff, "E", "D", true, "increase_ppi", false, 0.2);

		assertAnEdge(draughtDiff, "A", "B", false, "decreases_phosphorylation", false, 1.1);

		assertAnEdge(draughtDiff, "M", "N", false, "increases_phosphorylation", false, 4);
		assertAnEdge(draughtDiff, "P", "M", false, "increases_ptm", false, 7);
		assertAnEdge(draughtDiff, "N", "O", false, "increases_phosphorylation", false, 5);
		assertAnEdge(draughtDiff, "N", "P", false, "decreases_phosphorylation", false, 3);
		assertAnEdge(draughtDiff, "P", "N", false, "increases_phosphorylation", false, 8);

		// Testing the edges in the corresponding consensus network
		ConsensusNetwork draughtConsensus = draughtPair.getConsensusNetwork();
		Set<Edge> cEdges = draughtConsensus.getEdges();
		assertEquals(3, cEdges.size());

		assertAnEdge(draughtConsensus, "A", "B", true, "ppi", false, 0.3);
		assertAnEdge(draughtConsensus, "A", "C", true, "ppi", false, 0.8);

		assertAnEdge(draughtConsensus, "M", "N", false, "phosphorylation", false, 2);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results.
	 * This method specifically checks the pairwise consensus algorithms of a generic set of input networks (i.e. reference undefined).
	 */
	@Test
	public void testMultipleConditionsPairwiseConsensus()
	{
		MultipleConditionTest ex = new MultipleConditionTest();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int ID = ex.getDefaultRunConfigurationID(p);
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, weight_cutoff, false, true, 10, true, null);

		// Testing that there are exactly three consensus networks created (3 pairs)
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 0);
		assertNrConsensusNetworks(output, 3);

		Map<String, ConsensusNetwork> networks = new HashMap<String, ConsensusNetwork>();
		for (ConsensusNetwork on : output.getConsensusNetworks())
		{
			networks.put(on.getName(), on);
		}

		// Salt vs. reference
		ConsensusNetwork saltConsensus = networks.get("consensus_Reference_Salty");
		Set<Edge> sEdgesS = saltConsensus.getEdges();
		assertEquals(5, sEdgesS.size());

		assertAnEdge(saltConsensus, "A", "Z", true, "ppi", false, 0.1);
		assertAnEdge(saltConsensus, "A", "B", true, "ppi", false, 0.4);
		assertAnEdge(saltConsensus, "A", "C", true, "ppi", false, 0.6);

		assertAnEdge(saltConsensus, "M", "N", false, "phosphorylation", false, 2);
		assertAnEdge(saltConsensus, "M", "O", false, "phosphorylation", true, 1);

		// Draught vs. reference
		ConsensusNetwork draughtConsensus = networks.get("consensus_Draughty_Reference");
		Set<Edge> sEdges = draughtConsensus.getEdges();
		assertEquals(3, sEdges.size());

		assertAnEdge(draughtConsensus, "A", "B", true, "ppi", false, 0.3);
		assertAnEdge(draughtConsensus, "A", "C", true, "ppi", false, 0.8);

		assertAnEdge(draughtConsensus, "M", "N", false, "phosphorylation", false, 2);

		// Draught vs. reference
		ConsensusNetwork draughtStressConsensus = networks.get("consensus_Draughty_Salty");
		Set<Edge> sEdgesDS = draughtStressConsensus.getEdges();
		assertEquals(5, sEdgesDS.size());

		assertAnEdge(draughtStressConsensus, "A", "B", true, "ppi", false, 0.3);
		assertAnEdge(draughtStressConsensus, "A", "C", true, "ppi", false, 0.6);
		assertAnEdge(draughtStressConsensus, "A", "D", true, "ppi", false, 0.75);

		assertAnEdge(draughtStressConsensus, "M", "N", false, "phosphorylation", false, 6);
		assertAnEdge(draughtStressConsensus, "P", "M", false, "ptm", false, 2);
	}

	/**
	 * JUNIT Test: check whether the example network with edge conflicts produces correct results.
	 */
	@Test
	public void testConflicts()
	{
		ConflictingEdgesTest ex = new ConflictingEdgesTest();
		double cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int ID = ex.getDefaultRunConfigurationID(p);
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, cutoff, null, null, 10, 20, true, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 1);

		// Testing the edges in the differential network
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
		DifferentialNetwork dNetwork = pair.getDifferentialNetwork();

		Set<Edge> dEdges = dNetwork.getEdges();

		assertEquals(9, dEdges.size());

		assertAnEdge(dNetwork, "A", "B", false, "increases_regulation", false, 6);
		assertAnEdge(dNetwork, "A", "B", false, "decreases_ptm", false, 5);
		assertAnEdge(dNetwork, "A", "B", false, "decreases_somerandomInteraction", false, 4);

		assertAnEdge(dNetwork, "G", "H", false, "decreases_unspecified_regulation", false, 3);
		assertAnEdge(dNetwork, "G", "H", false, "decreases_ptm", false, 1);

		assertAnEdge(dNetwork, "J", "K", false, "decreases_unspecified_regulation", false, 3);
		assertAnEdge(dNetwork, "K", "J", false, "decreases_regulation", false, 4);
		assertAnEdge(dNetwork, "J", "K", false, "increases_ptm", false, 4);
		assertAnEdge(dNetwork, "K", "J", false, "increases_ptm", false, 1);

		// Testing the edges in the corresponding consensus network
		ConsensusNetwork sNetwork = pair.getConsensusNetwork();
		Set<Edge> sEdges = sNetwork.getEdges();
		assertEquals(7, sEdges.size());

		assertAnEdge(sNetwork, "A", "B", false, "positive_regulation", false, 2);

		assertAnEdge(sNetwork, "G", "H", false, "regulation", false, 4);
		assertAnEdge(sNetwork, "G", "H", false, "ptm", false, 2);

		assertAnEdge(sNetwork, "J", "K", false, "regulation", false, 4);
		assertAnEdge(sNetwork, "J", "K", false, "positive_regulation", true, 3);
		assertAnEdge(sNetwork, "K", "J", false, "ptm", false, 2);
		assertAnEdge(sNetwork, "J", "K", false, "ptm", false, 2);
	}
}
