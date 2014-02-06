package be.svlandeg.diffany.cytoscape.tasks;

import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskIterator;

import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.cytoscape.internal.Services;

/**
 * Factory to call {@link LoadExampleTask}
 * 
 * @author Thomas Van Parys
 *
 */
public class LoadExampleTaskFactory implements TaskFactory{

	private Services services;
	private Project exampleProject;
	
	/**
	 * 
	 * @param services the app {@link Services}
	 * @param exampleProject {@link Project} to be used as example input.
	 */
	public LoadExampleTaskFactory(Services services, Project exampleProject) {
		this.services = services;
		this.exampleProject = exampleProject;
	}

	@Override
	public TaskIterator createTaskIterator() {
		TaskIterator it = new TaskIterator();

		it.append(new LoadExampleTask(services, exampleProject));

		return it;
	}

	@Override
	public boolean isReady() {
		return true;
	}

}
