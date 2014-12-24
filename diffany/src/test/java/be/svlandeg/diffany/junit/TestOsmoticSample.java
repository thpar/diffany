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

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.examples.OsmoticSampleTest;

/** 
 * This class provides examples to benchmark the fuzzy consensus & differential functionality, using real examples from the osmotic data.
 * 
 * @author Sofie Van Landeghem
 */
public class TestOsmoticSample extends TestGeneric
{
	 
	/**
	 * JUNIT Test: check whether the example fuzzy network produces correct results with the min operator and cutoff 4 out of 4
	 */
	@Test
	public void testSample()
	{
		OsmoticSampleTest ex = new OsmoticSampleTest();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int supportingCutoff = 5;
		int ID = ex.getTestDiffConfiguration(p, supportingCutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, null, null, 70, 80, true, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 1);
		
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();

		// Testing the edges in the differential network
		DifferentialNetwork dn = pair.getDifferentialNetwork();
		Set<Edge> dEdges = dn.getEdges();
		assertEquals(2, dEdges.size());
		
		assertAnEdge(dn, "M", "N", false, "decreases_regulation", false, 1.1);
		assertAnEdge(dn, "O", "P", false, "decreases_regulation", false, 1.1);
	}

}
