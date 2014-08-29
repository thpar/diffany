package be.svlandeg.diffany.core.project;

import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.meta.MetaConvertor;
import be.svlandeg.diffany.core.networks.meta.MetaDifferentialNetwork;
import be.svlandeg.diffany.core.networks.meta.MetaConsensusNetwork;

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
			String errormsg = "The specified output differential network(s) can not be null!";
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
	 * @param on the consensus output network
	 */
	public void addConsensus(ConsensusNetwork cn)
	{
		if (cn == null)
		{
			String errormsg = "The specified output consensus network(s) can not be null!";
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

	/**
	 * Retrieve all differential networks as one large merged network
	 * @return the differential networks in this output, all merged together (or null when there were no differential networks)
	 */
	public MetaDifferentialNetwork getMergedDifferential()
	{
		if (dns == null || dns.isEmpty())
		{
			return null;
		}
		return MetaConvertor.convertDifferentials(dns);
	}

	/**
	 * Retrieve all consensus networks as one large merged network
	 * @return the consensus networks in this output, all merged together
	 */
	public MetaConsensusNetwork getMergedConsensusNetworks()
	{
		if (cns == null || cns.isEmpty())
		{
			return null;
		}
		return MetaConvertor.convertConsensus(cns);
	}

}
