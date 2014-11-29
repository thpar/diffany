package be.svlandeg.diffany.cytoscape.gui;

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
