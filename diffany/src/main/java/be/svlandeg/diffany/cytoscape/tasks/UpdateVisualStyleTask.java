package be.svlandeg.diffany.cytoscape.tasks;

import java.util.Set;

import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.vizmap.VisualStyle;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskMonitor;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.semantics.EdgeOntology;

/**
 * This task gathers all used interactions in a given {@link CyProject}, updates
 * this {@link VisualStyle}s and applies them to all {@link CyView}s in the {@link CyProject}.
 * 
 * @author thpar
 *
 */
public class UpdateVisualStyleTask implements Task {

	private CyProject project;
	private Model model;
	
	public UpdateVisualStyleTask(Model model, CyProject project) {		
		this.project = project;
		this.model = model;
	}

	@Override
	public void run(TaskMonitor taskMonitor) throws Exception {
		
		Set<String> interactions = project.getAllInteractions();
		EdgeOntology ontology = project.getEdgeOntology();
		
		model.getSourceStyle().updateInteractionMappings(interactions, ontology);
		model.getDiffStyle().updateInteractionMappings(interactions, ontology);
		
		for (CyNetworkView view : project.getAllSourceViews(model.getServices())){
			model.getSourceStyle().apply(view);
			view.updateView();
		}
		for (CyNetworkView view : project.getAllDifferentialViews(model.getServices())){
			model.getDiffStyle().apply(view);
			view.updateView();
		}
	}

	@Override
	public void cancel() {
		// TODO Auto-generated method stub

	}

}
