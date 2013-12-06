package be.svlandeg.diffany.cytoscape.tasks;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskMonitor;

import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.cytoscape.CyNetworkBridge;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.internal.Services;

/**
 * Temporary task to experiment.
 * 
 * @author thpar
 *
 */
public class TestTask implements Task {

	private Network network;
	private Model model;
	
	
	public TestTask(Network network, Model model) {
		this.network = network;
		this.model = model;
	}

	@Override
	public void cancel() {
		System.out.println("Task cancelled");
	}

	@Override
	public void run(TaskMonitor taskMonitor) throws Exception {
		taskMonitor.setTitle("Test Task");
		taskMonitor.setProgress(0.1);
		Thread.sleep(3000);
		taskMonitor.setProgress(0.6);
		Thread.sleep(3000);
		
		Services services = model.getServices();
		
		CyNetworkBridge bridge = new CyNetworkBridge(model);		
		CyNetwork cyNetwork = bridge.createCyNetwork(network);
		services.getCyNetworkManager().addNetwork(cyNetwork);
		CyNetworkView cyView = services.getCyNetworkViewFactory().createNetworkView(cyNetwork);
		services.getCyNetworkViewManager().addNetworkView(cyView);
		taskMonitor.setProgress(1.0);
	}

}
