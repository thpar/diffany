package be.svlandeg.diffany.cytoscape.gui;

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

import org.cytoscape.model.CyNetwork;

/**
 * Wrapper class for {@link CyNetwork}s to easily include them in GUI components by supplying a toString method returning 
 * the network name.
 * 
 * Also supplies extra options per network.
 * 
 * @author Thomas Van Parys
 *
 */
public class NetworkEntry {
	private CyNetwork network;
	private boolean isSelected;
	private boolean isReference;
	

	private static final String UNKNOWN = "network unknown";
	
	/**
	 * Returns a new entry, based on given {@link CyNetwork}
	 * @param network the given network
	 */
	public NetworkEntry(CyNetwork network) {
		this.network = network;
	}
	
	/**
	 * Gets the default name from the network table.
	 * 
	 * @return String default name from the network table.
	 */
	public String getName(){
		try{
			return network.getRow(network).get(CyNetwork.NAME, String.class);			
		} catch (NullPointerException e){
			System.err.println(e.getMessage());
			return UNKNOWN;
		}
	}
	@Override
	public String toString(){
		return getName();
	}
	/**
	 * The wrapped {@link CyNetwork}
	 * @return The wrapped {@link CyNetwork}
	 */
	public CyNetwork getNetwork() {
		return network;
	}
	/**
	 * 
	 * @return is the network selected to be included in the current project?
	 */
	public boolean isSelected() {
		return isSelected;
	}
	/**
	 * Toggle whether or not the network should be used with the next execution of the algorithm.
	 * @param isSelected whether or not this network is selected
	 */
	public void setSelected(boolean isSelected) {
		this.isSelected = isSelected;
	}
	/**
	 * 
	 * @return is this network the reference network in the current project?
	 */
	public boolean isReference() {
		return isReference;
	}
	/**
	 * Toggles this network as the reference network on the next execution of the algorithm.
	 * @param isReference whether or not this network is the reference network
	 */
	public void setReference(boolean isReference) {
		this.isReference = isReference;
	}
	
		
	
}
