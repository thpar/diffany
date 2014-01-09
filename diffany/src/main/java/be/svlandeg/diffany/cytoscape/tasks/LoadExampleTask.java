package be.svlandeg.diffany.cytoscape.tasks;

import java.util.Collection;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;
import org.cytoscape.view.layout.CyLayoutAlgorithm;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.CyNetworkViewFactory;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskManager;
import org.cytoscape.work.TaskMonitor;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.cytoscape.CyNetworkBridge;
import be.svlandeg.diffany.internal.Services;

public class LoadExampleTask implements Task{

	Services services;
	private Project exampleProject ;
	
	public LoadExampleTask(Services services, Project exampleProject) {
		this.services = services;
		this.exampleProject = exampleProject;
	}

	@Override
	public void run(TaskMonitor taskMonitor) throws Exception {
		
		ReferenceNetwork refNet = exampleProject.getReferenceNetwork();
		Collection<ConditionNetwork> condNets = exampleProject.getConditionNetworks();
		
		CyNetworkManager networkManager = services.getCyNetworkManager();
		CyRootNetworkManager rootNetworkManager = services.getCyRootNetworkManager();
		CyNetworkViewFactory viewFactory = services.getCyNetworkViewFactory();
		CyNetworkViewManager viewManager = services.getCyNetworkViewManager();
		CyLayoutAlgorithm layout = services.getCyLayoutAlgorithmManager().getLayout("force-directed");
		TaskManager tm = services.getTaskManager();
		CyNetworkBridge bridge = new CyNetworkBridge();
		
		CyNetwork cyRefNet = bridge.createCyNetwork(refNet, services.getCyNetworkFactory());
		networkManager.addNetwork(cyRefNet);
		CyRootNetwork collection = rootNetworkManager.getRootNetwork(cyRefNet);
		collection.getRow(collection).set(CyRootNetwork.NAME, exampleProject.getName());
		
		CyNetworkView refView = viewFactory.createNetworkView(cyRefNet);
		viewManager.addNetworkView(refView);
		TaskIterator it = layout.createTaskIterator(refView, layout.createLayoutContext(), 
				CyLayoutAlgorithm.ALL_NODE_VIEWS, null);
		tm.execute(it);
		
		for (ConditionNetwork condNet : condNets){
			CyNetwork cyCondNet = bridge.createCyNetwork(condNet, collection);
			networkManager.addNetwork(cyCondNet);
			CyNetworkView condView = viewFactory.createNetworkView(cyCondNet);
			viewManager.addNetworkView(condView);
			it = layout.createTaskIterator(condView, layout.createLayoutContext(), 
					CyLayoutAlgorithm.ALL_NODE_VIEWS, null);
			tm.execute(it);
		}
		
	}

	@Override
	public void cancel() {
		// TODO Auto-generated method stub
		
	}

}
