package be.svlandeg.diffany.cytoscape.tasks;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskMonitor;

import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.cytoscape.CyNetworkBridge;
import be.svlandeg.diffany.internal.Services;

public class TestTask implements Task {

	private Network network;
	private Services services;
	
	
	public TestTask(Network network, Services services) {
		this.network = network;
		this.services = services;
	}

	@Override
	public void cancel() {
		// TODO Auto-generated method stub

	}

	@Override
	public void run(TaskMonitor arg0) throws Exception {
		CyNetworkBridge bridge = new CyNetworkBridge(services.getCyNetworkFactory());
		CyNetwork cyNetwork = bridge.createCyNetwork(network);
		services.getCyNetworkManager().addNetwork(cyNetwork);
		CyNetworkView cyView = services.getCyNetworkViewFactory().createNetworkView(cyNetwork);
		services.getCyNetworkViewManager().addNetworkView(cyView);
	}

}
