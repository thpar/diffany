package be.svlandeg.diffany.cytoscape;

import java.util.HashSet;
import java.util.Observable;
import java.util.Set;

import javax.swing.JFrame;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.events.NetworkAddedEvent;
import org.cytoscape.model.events.NetworkAddedListener;
import org.cytoscape.model.events.NetworkDestroyedEvent;
import org.cytoscape.model.events.NetworkDestroyedListener;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;

import be.svlandeg.diffany.cytoscape.vizmapper.VisualDiffStyle;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualSourceStyle;
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
public class Model extends Observable implements NetworkAddedListener, NetworkDestroyedListener{

	/**
	 * A collection of all Cytoscape services that were registered in the {@link CyActivator}
	 */
	private Services services;
	
	/**
	 * The current project
	 */
	private CyProject currentProject = new CyProject();

	/**
	 * Network collection that's currently selected
	 */
	private CyRootNetwork selectedCollection;
	
	private VisualSourceStyle sourceStyle;
	private VisualDiffStyle diffStyle;

	private JFrame swingApplication;
	
	
	public Model(Services services){
		this.services = services;
	
		sourceStyle = new VisualSourceStyle(services);
		diffStyle = new VisualDiffStyle(services);
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
		
		//display set
		System.out.println("Found network collections:");
		for (CyRootNetwork net: set){
			String name = net.getRow(net).get(CyNetwork.NAME, String.class);
			System.out.println(" - "+name);
		}
		return set;
	}
	
	/**
	 * Returns the {@link CyNetwork} with given name.
	 * Returns null if no such network exists.
	 * 
	 * @param id the network name
	 * @return {@link CyNetwork} with given id
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
		setChanged();
		notifyObservers();
	}


	@Override
	public void handleEvent(NetworkAddedEvent e) {
		//triggered on network added
		setChanged();
		notifyObservers();
	}
	
	@Override
	public void handleEvent(NetworkDestroyedEvent e) {
		currentProject.removeDestroyedNetworks();
		setChanged();
		notifyObservers();
	}
	

	public VisualSourceStyle getSourceStyle() {
		return sourceStyle;
	}


	public VisualDiffStyle getDiffStyle() {
		return diffStyle;
	}


	public void setParentWindow(JFrame jFrame) {
		this.swingApplication = jFrame;
	}

	public JFrame getParentWindow(){
		return this.swingApplication;
	}

	
}
