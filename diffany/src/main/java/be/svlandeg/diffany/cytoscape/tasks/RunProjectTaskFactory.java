package be.svlandeg.diffany.cytoscape.tasks;

import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskIterator;

import be.svlandeg.diffany.cytoscape.Model;

/**
 * Temporary {@link TaskFactory} to execute the {@link TestTask}.
 * 
 * @author thpar
 *
 */
public class RunProjectTaskFactory implements TaskFactory {

	private Model model;
	
	public RunProjectTaskFactory(Model model) {
		this.model = model;
	}

	@Override
	public TaskIterator createTaskIterator() {
		TaskIterator it = new TaskIterator(new RunProjectTask(model));
		return it;
	}

	@Override
	public boolean isReady() {
		return model.getCurrentProject().canExecute();
	}

}
