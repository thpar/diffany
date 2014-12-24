package be.svlandeg.diffany.core.project;

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


import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;

/**
 * This class keeps the output of differential network algorithms, both differential networks as well consensus networks.
 * Either of the two can be null when it was not calculated.
 * 
 * @author Sofie Van Landeghem
 */
public class RunOutput
{
	
	protected Project p;
	protected int runID;

	private Set<DifferentialNetwork> dns;
	private Set<ConsensusNetwork> cns;
	private Set<OutputNetworkPair> pairs;

	/**
	 * Create a new empty output object.
	 * @param p the project to which this output object belongs
	 * @param runID the corresponding run ID within the project
	 *
	 * @throws IllegalArgumentException when either of the provided networks is null
	 */
	public RunOutput(Project p, int runID)
	{
		clean();
		this.p = p;
		this.runID = runID;
	}

	/**
	 * Clean the output.
	 */
	public void clean()
	{
		dns = new HashSet<DifferentialNetwork>();
		cns = new HashSet<ConsensusNetwork>();
		pairs = new HashSet<OutputNetworkPair>();
	}

	/**
	 * Add a pair of differential+consensus networks to this output. 
	 * This will automatically also register the output networks to the list of differential and consensus networks, correspondingly.
	 * When adding these networks, the IDs of the output networks will be checked for uniquenesss.
	 * 
	 * @param pair the pair of differential and consensus output networks
	 */
	public void addPair(OutputNetworkPair pair)
	{
		if (pair == null)
		{
			String errormsg = "The specified output pair can not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		
		addDifferential(pair.getDifferentialNetwork());
		addConsensus(pair.getConsensusNetwork());
		
		pairs.add(pair);
	}

	/**
	 * Add a differential output network. When adding this network, its IDs will be checked for uniquenesss.
	 * @param dn the differential output network
	 */
	public void addDifferential(DifferentialNetwork dn)
	{
		if (dn == null)
		{
			String errormsg = "The specified output differential network can not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		
		if (! p.runs.get(runID).checkoutputID(dn.getID()))
		{
			String errormsg = "The provided ID of the differential network (" + dn.getID() + ") should be unique in this run!";
			throw new IllegalArgumentException(errormsg);
		}
		
		dns.add(dn);
	}

	/**
	 * Add a consensus output network. When adding this network, its IDs will be checked for uniquenesss.
	 * @param cn the consensus output network
	 */
	public void addConsensus(ConsensusNetwork cn)
	{
		if (cn == null)
		{
			String errormsg = "The specified output consensus network can not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		
		if (! p.runs.get(runID).checkoutputID(cn.getID()))
		{
			String errormsg = "The provided ID of the consensus network (" + cn.getID() + ") should be unique in this run!";
			throw new IllegalArgumentException(errormsg);
		}
		
		cns.add(cn);
	}

	/**
	 * Retrieve the differential and consensus networks as result pairs.
	 * @return the result pairs, each containing both a differential and a consensus network
	 */
	public Set<OutputNetworkPair> getOutputAsPairs()
	{
		return pairs;
	}

	/**
	 * Retrieve all differential networks
	 * @return all differential networks in this output
	 */
	public Set<DifferentialNetwork> getDifferentialNetworks()
	{
		return dns;
	}

	/**
	 * Retrieve all consensus networks
	 * @return the consensus network in this output
	 */
	public Set<ConsensusNetwork> getConsensusNetworks()
	{
		return cns;
	}

}
