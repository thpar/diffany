package be.svlandeg.diffany.cytoscape.tasks;

import java.util.Set;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.vizmap.VisualStyle;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskMonitor;

import be.svlandeg.diffany.core.semantics.EdgeOntology;
import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;

/**
 * This task gathers all used interactions in a given {@link CyProject}, updates
 * this {@link VisualStyle}s and applies them to all {@link CyNetworkView}s in the {@link CyProject}.
 * 
 * @author Thomas Van Parys
 *
 */
public class UpdateVisualStyleTask implements Task {

	private CyProject cyProject;
	private Model model;
	
	/**
	 * Construct a new task to update the visual styles in the model, according to the interactions 
	 * in the {@link CyNetwork}s contained in the {@link CyProject}. The new visual styles will be applied 
	 * to the same set of networks.
	 * 
	 * @param model Diffany {@link Model}
	 * @param project the {@link CyProject} containing the used {@link CyNetwork}s. 
	 */
	public UpdateVisualStyleTask(Model model, CyProject project) {		
		this.cyProject = project;
		this.model = model;
	}

	@Override
	public void run(TaskMonitor taskMonitor) throws Exception {
		
		Set<String> sourceInteractions = cyProject.getAllSourceInteractions();
		Set<String> diffInteractions = cyProject.getAllDifferentialInteractions();
		EdgeOntology ontology = cyProject.getProject().getEdgeOntology();
		
		model.getSourceStyle().updateInteractionMappings(sourceInteractions, ontology);
		model.getDiffStyle().updateInteractionMappings(diffInteractions, ontology);
		
		for (CyNetworkView view : cyProject.getAllSourceViews(model.getServices())){
			model.getSourceStyle().apply(view);
			view.updateView();
		}
		for (CyNetworkView view : cyProject.getAllDifferentialViews(model.getServices())){
			model.getDiffStyle().apply(view);
			view.updateView();
		}
	}

	@Override
	public void cancel() {
		// TODO Auto-generated method stub

	}

}
