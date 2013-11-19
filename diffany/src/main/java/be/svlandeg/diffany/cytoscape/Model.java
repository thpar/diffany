package be.svlandeg.diffany.cytoscape;

import java.util.Observable;

import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Project;
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
	
	public Model(Services services){
		this.services = services;
	}
	
	/** 
	 * Use information from this model to construct the {@link Project} and 
	 * {@link Network}s to run the actual algorithm. 
	 */
	public void runAlgorithm(){
		
	}

	public Services getServices() {
		return services;
	}
}
