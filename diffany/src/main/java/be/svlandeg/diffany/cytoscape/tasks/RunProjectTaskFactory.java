package be.svlandeg.diffany.cytoscape.tasks;

import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskIterator;

import be.svlandeg.diffany.cytoscape.Model;

/**
 * Factory to run {@link RunProjectTask}, followed by a {@link UpdateVisualStyleTask}
 * 
 * @author Thomas Van Parys
 *
 */
public class RunProjectTaskFactory implements TaskFactory {

	private Model model;
	
	/**
	 * Construct a new factory.
	 * 
	 * @param model Diffany {@link Model}
	 */
	public RunProjectTaskFactory(Model model) {
		this.model = model;
	}

	@Override
	public TaskIterator createTaskIterator() {
		TaskIterator it = new TaskIterator();
		
		it.append(new RunProjectTask(model, model.getSelectedProject()));
		it.append(new UpdateVisualStyleTask(model, model.getSelectedProject()));
		
		return it;
	}

	@Override
	public boolean isReady() {
		return model.getSelectedProject().canExecute(model);
	}

}
