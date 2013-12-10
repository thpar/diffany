package be.svlandeg.diffany.cytoscape.vizmapper;

import java.util.HashSet;
import java.util.Set;

import org.cytoscape.view.vizmap.VisualStyle;

import be.svlandeg.diffany.internal.Services;

public abstract class VisualDiffanyStyle {

	protected String name;
	protected VisualStyle vis;
	protected Services services;
	
	/**
	 * Create a new visual style and initialize it according to its type. Then register the
	 * new style with Cytoscape.
	 * 
	 * @param name the displayed name in the VizMapper
	 * @param services the services object
	 */
	public VisualDiffanyStyle(String name, Services services){
		this.name = name;
		this.services = services;
		this.vis = services.getVisualStyleFactory().createVisualStyle(name);
		services.getVisualMappingManager().addVisualStyle(this.vis);
		
		this.init();
	}

	
	abstract protected void init();
	
	public VisualStyle getVisualStyle(){
		return this.vis;
	}
	
	private Set<String> getAllInteractions(){
		Set<String> interactions = new HashSet<String>();
		
		return interactions;
	}
}
