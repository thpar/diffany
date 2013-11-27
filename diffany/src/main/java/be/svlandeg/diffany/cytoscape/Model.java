package be.svlandeg.diffany.cytoscape;

import java.util.HashSet;
import java.util.Observable;
import java.util.Set;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;

import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.cytoscape.gui.GUIModel;
import be.svlandeg.diffany.internal.Services;

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
	
	public Model(Services services){
		this.services = services;
		this.guiModel = new GUIModel();
	}
	
	public GUIModel getGuiModel(){
		return this.guiModel;
	}
	
	/** 
	 * Use information from this model to construct the {@link Project} and 
	 * {@link Network}s to run the actual algorithm. 
	 */
	public void runAlgorithm(){
		
	}
	
		

	public Set<CyRootNetwork> getNetworkCollections(){
		Set<CyNetwork> allNetworks = services.getCyNetworkManager().getNetworkSet();
		CyRootNetworkManager rootNetManager = services.getCyRootNetworkManager();
		
		Set<CyRootNetwork> set = new HashSet<CyRootNetwork>();
		
		for (CyNetwork net : allNetworks){
			set.add(rootNetManager.getRootNetwork(net));
		}
		
		return set;
	}
	
	public Services getServices() {
		return services;
	}
}
