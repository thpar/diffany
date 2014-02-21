package be.svlandeg.diffany.cytoscape.tasks;

import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskIterator;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;

/**
 * Factory to run the {@link UpdateVisualStyleTask}
 * 
 * @author Thomas Van Parys
 *
 */
public class UpdateVisualStyleTaskFactory implements TaskFactory {

	private Model model;
	private CyProject cyProject;
	
	/**
	 * Construct a new factory to create {@link UpdateVisualStyleTask}s
	 * 
	 * @param model Diffany {@link Model}
	 * @param project the {@link CyProject} containing the used {@link CyNetwork}s.
	 */
	public UpdateVisualStyleTaskFactory(Model model, CyProject project) {
		this.model = model;
		this.cyProject = project;
	}

	@Override
	public TaskIterator createTaskIterator() {
		TaskIterator it = new TaskIterator(new UpdateVisualStyleTask(model, cyProject));
		return it;
	}

	@Override
	public boolean isReady() {
		return true;
	}

}
