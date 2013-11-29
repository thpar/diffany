package be.svlandeg.diffany.cytoscape.gui;

import java.util.ArrayList;
import java.util.List;
import java.util.Observable;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CySubNetwork;

import be.svlandeg.diffany.cytoscape.NetworkEntry;

public class GUIModel extends Observable{

	private CyNetwork selectedCollection;
	
	private CyNetwork referenceNetwork;
	private List<NetworkEntry> networkEntries = new ArrayList<NetworkEntry>();
	
		
	public CyNetwork getSelectedCollection() {
		return selectedCollection;
	}

	public void setSelectedCollection(CyNetwork selectedCollection) {
		this.selectedCollection = selectedCollection;
		refreshNetworkEntries();
		setChanged();
		notifyObservers("Boom!!!");
	}

	public List<NetworkEntry> getNetworkEntries(){
		return this.networkEntries;
	}
	
	private void refreshNetworkEntries() {
		List<CySubNetwork> subNets = ((CyRootNetwork)selectedCollection).getSubNetworkList();
		networkEntries = new ArrayList<NetworkEntry>();
		for (CySubNetwork subNet : subNets){
			NetworkEntry entry = new NetworkEntry();
			entry.setNetwork(subNet);
			entry.setSelected(false);
			entry.setReference(false);
			networkEntries.add(entry);
		}
		if (networkEntries.size() > 0){
			networkEntries.get(0).setReference(true);
			this.referenceNetwork = networkEntries.get(0).getNetwork();
		}
	}

	public CyNetwork getReferenceNetwork() {
		return referenceNetwork;
	}
	
		
}
