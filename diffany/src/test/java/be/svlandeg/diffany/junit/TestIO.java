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
import static org.junit.Assert.fail;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.io.NetworkIO;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunDiffConfiguration;
import be.svlandeg.diffany.examples.Bandyopadhyay2010;

/** 
 * Class that automatically tests the IO functionality of Diffany:
 * writing/reading networks and projects to/from file.
 * 
 * @author Sofie Van Landeghem
 */
public class TestIO
{
	
	/**
	 * JUNIT Test: write and read all types of Networks.
	 * The networks used in this tests are tested for completeness in {@link TestExamples#testBandyopadhyay}
	 */
	@Test
	public void testNetworkIO()
	{
		// System-dependent tmp dir. E.g. windows 7: C:\Users\YourUserName\AppData\Local\Temp
		String testLocation = System.getProperty("java.io.tmpdir") + File.separator + "diffany" + File.separator;
		File rDir = new File(testLocation + "reference/");
		File cDir = new File(testLocation + "condition/");
		File dDir = new File(testLocation + "differential/");
		File consDir = new File(testLocation + "consensus/");
		
		rDir.mkdirs();
		cDir.mkdirs();
		dDir.mkdirs();
		consDir.mkdirs();
		
		Bandyopadhyay2010 ex = new Bandyopadhyay2010();
		Project p = ex.getDefaultProject();
		int ID = ex.getDefaultRunConfigurationID(p);
		double cutoff = 0.0;
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, cutoff, true, true, 20, true, null);
		
		RunDiffConfiguration rc = (RunDiffConfiguration) p.getRunConfiguration(ID);
		
		// The reference network
		ReferenceNetwork rWriteNetwork = rc.getReferenceNetwork();
		
		// The condition-dependent network (there is only 1 in that specific project)
		ConditionNetwork cWriteNetwork = rc.getConditionNetworks().iterator().next();
		
		// There is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
		DifferentialNetwork dNetwork = pair.getDifferentialNetwork();
				
		// The consensus network (there should be only 1)
		ConsensusNetwork cNetwork = pair.getConsensusNetwork();
			
		// WRITING
		boolean writeHeader = true;
		try
		{
			NetworkIO.writeNetworkToDir(rWriteNetwork, rDir, writeHeader);
			NetworkIO.writeNetworkToDir(cWriteNetwork, cDir, writeHeader);
			NetworkIO.writeNetworkToDir(dNetwork, dDir, writeHeader);
			NetworkIO.writeNetworkToDir(cNetwork, consDir, writeHeader);
		}
		catch(IOException io)
		{
			fail(io.getMessage());
		}
		
		// READING
		try
		{
			ReferenceNetwork rReadNetwork = NetworkIO.readReferenceNetworkFromDir(rDir, writeHeader, true);
			assertEquals(3, rReadNetwork.getEdges().size());
			assertEquals(5, rReadNetwork.getNodes().size());
			
			ConditionNetwork cReadNetwork = NetworkIO.readConditionNetworkFromDir(cDir, writeHeader, true);
			assertEquals(3, cReadNetwork.getEdges().size());
			assertEquals(4, cReadNetwork.getNodes().size());
			
			Set<ConditionNetwork> cReadNetworks = new HashSet<ConditionNetwork>();
			cReadNetworks.add(cReadNetwork);
			
			DifferentialNetwork dReadNetwork = NetworkIO.readDifferentialNetworkFromDir(dDir, rReadNetwork, cReadNetworks, writeHeader);
			assertEquals(3, dReadNetwork.getEdges().size());
			assertEquals(4, dReadNetwork.getNodes().size());
			
			ConsensusNetwork consReadNetwork = NetworkIO.readConsensusNetworkFromDir(consDir, rReadNetwork, cReadNetworks, writeHeader);
			assertEquals(2, consReadNetwork.getEdges().size());
			assertEquals(3, consReadNetwork.getNodes().size());
		}
		catch(IOException io)
		{
			fail(io.getMessage());
		}
		
	}
}
