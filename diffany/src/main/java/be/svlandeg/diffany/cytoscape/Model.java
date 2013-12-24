package be.svlandeg.diffany.cytoscape;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Observable;
import java.util.Set;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;
import org.cytoscape.model.subnetwork.CySubNetwork;

import be.svlandeg.diffany.internal.CyActivator;
import be.svlandeg.diffany.internal.Services;

/**
 * Model that keeps track of all settings and selections within the Cytoscape App.
 * Only calling {@link runAlgorithm} should construct the actual model to do the 
 * calculations and produce results, which are handed back to this model.
 * 
 * @author thpar
 *
 */
public class Model extends Observable {

	/**
	 * A collection of all Cytoscape services that were registered in the {@link CyActivator}
	 */
	private Services services;
	
	/**
	 * The current project
	 */
	private CyProject currentProject = new CyProject(this);

	/**
	 * Network collection that's currently selected
	 */
	private CyRootNetwork selectedCollection;
	
	/**
	 * Networks to be listed in the selection list
	 */
	private List<NetworkEntry> networkEntries = new ArrayList<NetworkEntry>();
	
	
	public Model(Services services){
		this.services = services;
	
	}
	

	/**
	 * Iterates over all networks in the current cytoscape session and returns the set of network collections.
	 * 
	 * @return the set of network collections.
	 */
	public Set<CyRootNetwork> getNetworkCollections(){
		Set<CyNetwork> allNetworks = services.getCyNetworkManager().getNetworkSet();
		CyRootNetworkManager rootNetManager = services.getCyRootNetworkManager();
		
		Set<CyRootNetwork> set = new HashSet<CyRootNetwork>();
		
		for (CyNetwork net : allNetworks){
			set.add(rootNetManager.getRootNetwork(net));
		}
		
		return set;
	}
	
	/**
	 * Returns the {@link CyNetwork} with given name.
	 * Returns null if no such network exists.
	 * 
	 * @param id the network name
	 * @return 
	 */
	public CyNetwork getNetworkByName(String id){
		Set<CyNetwork> allNetworks = services.getCyNetworkManager().getNetworkSet();
		for (CyNetwork net : allNetworks){
			String name = net.getRow(net).get(CyNetwork.NAME, String.class);
			if (name.equals(id)){
				return net;
			}
		}
		return null;
	}
	
	/**
	 * Gives access to a subset of the services offered in this context, as loaded in the {@link CyActivator}
	 * 
	 * @return
	 */
	public Services getServices() {
		return services;
	}

	
	/**
	 * Returns the current {@link CyProject}.
	 * @return the current {@link CyProject}
	 */
	public CyProject getCurrentProject() {
		return currentProject;
	}
	
	/**
	 * Erases all project settings and creates a fresh one.
	 * @return the newly created project
	 */
	public CyProject resetProject(){
		this.currentProject = new CyProject(this);
		return currentProject;
	}

	
	
	
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
		CyProject newProject = this.resetProject();
		
		
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


	
}
