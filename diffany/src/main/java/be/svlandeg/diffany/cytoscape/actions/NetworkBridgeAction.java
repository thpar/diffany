package be.svlandeg.diffany.cytoscape.actions;

import java.awt.event.ActionEvent;

import org.cytoscape.application.swing.AbstractCyAction;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskManager;
import org.cytoscape.work.swing.DialogTaskManager;

import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.tasks.TestTaskFactory;
import be.svlandeg.diffany.internal.Services;

public class NetworkBridgeAction extends AbstractCyAction {

	private Network network;
	private Services services;

	public NetworkBridgeAction(Model model, String menuTitle, Network network) {
		super(menuTitle, model.getServices().getCyApplicationManager(), null, null);
		setPreferredMenu("Apps.Diffany");
		this.services = model.getServices();
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

			TaskIterator it = tf.createTaskIterator();			
			DialogTaskManager dtm = services.getDialogTaskManager();
			dtm.execute(it);
		}
	}

}
