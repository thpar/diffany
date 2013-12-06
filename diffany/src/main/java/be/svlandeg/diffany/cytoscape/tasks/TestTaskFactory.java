package be.svlandeg.diffany.cytoscape.tasks;

import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskIterator;

import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.cytoscape.Model;

/**
 * Temporary {@link TaskFactory} to execute the {@link TestTask}.
 * 
 * @author thpar
 *
 */
public class TestTaskFactory implements TaskFactory {

	private Network network;
	private Model model;
	
	public TestTaskFactory(Network network, Model model) {
		this.network = network;
		this.model = model;
	}

	@Override
	public TaskIterator createTaskIterator() {
		TaskIterator it = new TaskIterator(new TestTask(network, model));
		return it;
	}

	@Override
	public boolean isReady() {
		return true;
	}

}
