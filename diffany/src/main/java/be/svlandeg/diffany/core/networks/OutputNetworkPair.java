package be.svlandeg.diffany.core.networks;

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



/**
 * This class simply contains one differential network and its corresponding consensus network.
 * 
 * @author Sofie Van Landeghem
 */
public class OutputNetworkPair
{
	
	private DifferentialNetwork dn;
	private ConsensusNetwork cn;
	
	/**
	 * Create a new pair from one specific differential network and its counterpart consensus network. Neither of the two networks may be null!
     *
	 * @param dn the differential network
	 * @param cn the consensus network
	 * @throws IllegalArgumentException when either of the provided networks is null
	 */
	public OutputNetworkPair(DifferentialNetwork dn, ConsensusNetwork cn)
	{
		if (dn == null || cn == null)
		{
			String errormsg = "The networks in the output network pair can not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.dn = dn;
		this.cn = cn;
	}
	
	/**
	 * Retrieve the differential network
	 * @return the differential network
	 */
	public DifferentialNetwork getDifferentialNetwork()
	{
		return dn;
	}
	
	/**
	 * Retrieve the consensus network
	 * @return the consensus network
	 */
	public ConsensusNetwork getConsensusNetwork()
	{
		return cn;
	}

}
