package be.svlandeg.diffany.cytoscape;

import java.util.Collection;
import java.util.HashSet;
import java.util.Observable;
import java.util.Set;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;
import org.cytoscape.view.model.CyNetworkView;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.cytoscape.gui.GUIModel;
import be.svlandeg.diffany.internal.CyActivator;
import be.svlandeg.diffany.internal.Services;
import be.svlandeg.diffany.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.semantics.DefaultNodeMapper;

/**
 * Model that keeps track of all settings and selections within the Cytoscape App.
 * Only calling {@link runAlgorithm} should construct the actual model to do the 
 * calculations and produce results, which are handed back to this model.
 * 
 * @author thpar
 *
 */
public class Model extends Observable{

	Services services;
	private GUIModel guiModel;
	
	private Project currentProject;
	
	public Model(Services services){
		this.services = services;
		this.guiModel = new GUIModel();
	}
	
	/**
	 * 
	 * @return the guiModel for this app.
	 */
	public GUIModel getGuiModel(){
		return this.guiModel;
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

	
	public Project getCurrentProject(){
		return currentProject;
	}
	
	
	
}
