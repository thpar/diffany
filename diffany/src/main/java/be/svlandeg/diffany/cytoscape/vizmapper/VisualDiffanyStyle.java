package be.svlandeg.diffany.cytoscape.vizmapper;

import org.cytoscape.view.vizmap.VisualStyle;

import be.svlandeg.diffany.internal.Services;

public abstract class VisualDiffanyStyle {

	protected String name;
	protected VisualStyle vis;
	
	/**
	 * Create a new visual style and initialize it according to its type. Then register the
	 * new style with Cytoscape.
	 * 
	 * @param name the displayed name in the VizMapper
	 * @param services the services object
	 */
	public VisualDiffanyStyle(String name, Services services){
		this.name = name;
		this.vis = services.getVisualStyleFactory().createVisualStyle(name);
		services.getVisualMappingManager().addVisualStyle(this.vis);
		this.init();
	}

	
	abstract protected void init();
}
