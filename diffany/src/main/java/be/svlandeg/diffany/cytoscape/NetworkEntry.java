package be.svlandeg.diffany.cytoscape;

import org.cytoscape.model.CyNetwork;

public class NetworkEntry {
	private CyNetwork network;
	private boolean isSelected;
	private boolean isReference;
	
	public NetworkEntry() {
	}
	public NetworkEntry(CyNetwork network) {
		this.network = network;
	}
	public String getName() {
		return network.getRow(network).get(CyNetwork.NAME, String.class);
	}
	@Override
	public String toString(){
		return getName();
	}
	public CyNetwork getNetwork() {
		return network;
	}
	public void setNetwork(CyNetwork network) {
		this.network = network;
	}
	public boolean isSelected() {
		return isSelected;
	}
	public void setSelected(boolean isSelected) {
		this.isSelected = isSelected;
	}
	public boolean isReference() {
		return isReference;
	}
	public void setReference(boolean isReference) {
		this.isReference = isReference;
	}
	
	
	
}
