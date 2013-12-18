package be.svlandeg.diffany.cytoscape.gui;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Observable;
import java.util.Set;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CySubNetwork;
import org.cytoscape.view.vizmap.VisualStyle;

import be.svlandeg.diffany.cytoscape.NetworkEntry;

/**
 * The GUIModel keeps track of the status of the different side panel components.
 * 
 * @author thpar
 *
 */
public class GUIModel extends Observable{

	private CyRootNetwork selectedCollection;
	private NetworkEntry referenceEntry;
	
	private List<NetworkEntry> networkEntries = new ArrayList<NetworkEntry>();
	
	private VisualStyle sourceStyle;
	private VisualStyle diffStyle;
	
	/**
	 * The collection of networks (aka the {@link CyRootNetwork}) that will be
	 * used for the algorithm.
	 * @return
	 */
	public CyRootNetwork getSelectedCollection() {
		return selectedCollection;
	}

	/**
	 * Set the collection of networks (aka the {@link CyRootNetwork}) that will be
	 * used for the algorithm.
	 * 
	 * This change will trigger an update with all observers.
	 * 
	 * @param selectedCollection
	 */
	public void setSelectedCollection(CyRootNetwork selectedCollection) {
		this.selectedCollection = selectedCollection;
		refreshNetworkEntries();
		setChanged();
		notifyObservers("Boom!!!");
	}

	/**
	 * The {@link NetworkEntry}s currently shown in the network selection table. These entries are wrappers for
	 * the list of {@link CySubNetwork}s from the selected network collection.
	 * @return current list entries in the network selection list
	 */
	public List<NetworkEntry> getNetworkEntries(){
		return this.networkEntries;
	}
	
	/**
	 * Reload the {@link NetworkEntry}s based on the selected network collection.
	 */
	private void refreshNetworkEntries() {
		List<CySubNetwork> subNets = ((CyRootNetwork)selectedCollection).getSubNetworkList();
		networkEntries = new ArrayList<NetworkEntry>();
		for (CySubNetwork subNet : subNets){
			NetworkEntry entry = new NetworkEntry(subNet);
			entry.setSelected(true);
			entry.setReference(false);
			networkEntries.add(entry);
		}
		if (networkEntries.size() > 0){
			networkEntries.get(0).setReference(true);
			this.referenceEntry = networkEntries.get(0);
		}
	}

	/**
	 * The currently selected reference network.
	 * @return
	 */
	public CyNetwork getReferenceNetwork() {
		return referenceEntry.getNetwork();
	}
	
	public Set<CyNetwork> getConditionEntries(){
		Set<CyNetwork> conditionals = new HashSet<CyNetwork>();
		for (NetworkEntry entry : this.networkEntries){
			if (entry.isSelected() && !entry.isReference()){
				conditionals.add(entry.getNetwork());
			}
		}
		return conditionals;
	}
	
	public NetworkEntry getReferenceEntry(){
		return referenceEntry;
	}
	
	public void setReferenceEntry(NetworkEntry entry){
		this.referenceEntry.setReference(false);
		this.referenceEntry = entry;
		entry.setReference(true);
	}

	public VisualStyle getSourceStyle() {
		return sourceStyle;
	}

	public void setSourceStyle(VisualStyle sourceStyle) {
		this.sourceStyle = sourceStyle;
	}

	public VisualStyle getDiffStyle() {
		return diffStyle;
	}

	public void setDiffStyle(VisualStyle diffStyle) {
		this.diffStyle = diffStyle;
	}
	
	
		
}
