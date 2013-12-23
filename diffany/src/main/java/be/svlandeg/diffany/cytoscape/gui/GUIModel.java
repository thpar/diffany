package be.svlandeg.diffany.cytoscape.gui;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Observable;
import java.util.Set;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CySubNetwork;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.NetworkEntry;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualDiffStyle;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualSourceStyle;

/**
 * The GUIModel keeps track of the status of the different side panel components.
 * 
 * @author thpar
 *
 */
public class GUIModel extends Observable{

	private Model model;
	
	public GUIModel(Model model){
		this.model = model;
	}
	
	
	private CyRootNetwork selectedCollection;
	
	/**
	 * Networks to be listed in the selection list
	 */
	private List<NetworkEntry> networkEntries = new ArrayList<NetworkEntry>();
	
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
		notifyObservers();
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
		
		//for now, we simply create a whole new CyProject
		//should become more flexible
		CyProject newProject = model.resetProject();
		
		
		for (CySubNetwork subNet : subNets){
			NetworkEntry entry = new NetworkEntry(subNet);
			entry.setSelected(false);
			entry.setReference(false);
			networkEntries.add(entry);
		}
		if (networkEntries.size() > 0){
			//set the table
			networkEntries.get(0).setReference(true);
			newProject.setReferenceNetwork(networkEntries.get(0).getNetwork());
			
		}
	}


	public Model getModel() {
		return model;
	}
		
		
}
