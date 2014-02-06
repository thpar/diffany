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
	private int runConfigurationID;
	
	/**
	 * 
	 * @param services the app {@link Services}
	 * @param exampleProject {@link Project} to be used as example input.
	 * @param runConfigurationID id of the example configuration to run in the project
	 */
	public LoadExampleTaskFactory(Services services, Project exampleProject, int runConfigurationID) {
		this.services = services;
		this.exampleProject = exampleProject;
		this.runConfigurationID = runConfigurationID;
	}

	@Override
	public TaskIterator createTaskIterator() {
		TaskIterator it = new TaskIterator();

		it.append(new LoadExampleTask(services, exampleProject, runConfigurationID));

		return it;
	}

	@Override
	public boolean isReady() {
		return true;
	}

}
