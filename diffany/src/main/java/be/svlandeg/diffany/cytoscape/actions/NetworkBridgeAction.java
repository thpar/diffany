package be.svlandeg.diffany.cytoscape.actions;

import java.awt.event.ActionEvent;

import org.cytoscape.application.swing.AbstractCyAction;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskManager;

import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.cytoscape.tasks.TestTaskFactory;
import be.svlandeg.diffany.internal.Services;

public class NetworkBridgeAction extends AbstractCyAction {

	private Network network;
	private Services services;

	public NetworkBridgeAction(Services services, String menuTitle, Network network) {
		super(menuTitle, services.getCyApplicationManager(), null, null);
		setPreferredMenu("Apps.Diffany");
		this.services = services;
		this.network = network;
	}

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	@Override
	public void actionPerformed(ActionEvent e) {
		TestTaskFactory tf = new TestTaskFactory(network, services);
		if (tf.isReady()){
			TaskManager<?, ?> tm = services.getTaskManager();
			TaskIterator it = tf.createTaskIterator();
			tm.execute(it);
		}
	}

}
