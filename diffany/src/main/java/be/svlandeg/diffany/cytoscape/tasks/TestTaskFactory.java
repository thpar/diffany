package be.svlandeg.diffany.cytoscape.tasks;

import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskIterator;

import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.internal.Services;

/**
 * Temporary {@link TaskFactory} to execute the {@link TestTask}.
 * 
 * @author thpar
 *
 */
public class TestTaskFactory implements TaskFactory {

	private Network network;
	private Services services;
	
	public TestTaskFactory(Network network, Services services) {
		this.network = network;
		this.services = services;
	}

	@Override
	public TaskIterator createTaskIterator() {
		TaskIterator it = new TaskIterator(new TestTask(network, services));
		return it;
	}

	@Override
	public boolean isReady() {
		return true;
	}

}
