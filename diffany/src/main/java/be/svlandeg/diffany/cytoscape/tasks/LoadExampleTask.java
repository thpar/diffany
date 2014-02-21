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
import be.svlandeg.diffany.cytoscape.internal.Services;

/**
 * Task to load example {@link Project}s into Cytoscape.
 * 
 * @author Thomas Van Parys
 *
 */
public class LoadExampleTask implements Task{

	private Services services;
	private Project exampleProject ;
	private int runConfigurationID;
	
	/**
	 * Construct a new task from to load networks from a {@link Project}. The resulting networks and 
	 * project settings will be ignored. Only the source networks are loaded and constructed. 
	 * 
	 * @param services the app {@link Services}
	 * @param exampleProject {@link Project} to be used as example input.
	 * @param runConfigurationID id of the example configuration to run in the project
	 */
	public LoadExampleTask(Services services, Project exampleProject, int runConfigurationID) {
		this.services = services;
		this.exampleProject = exampleProject;
		this.runConfigurationID = runConfigurationID;
	}

	@Override
	public void run(TaskMonitor taskMonitor) throws Exception {
		
		ReferenceNetwork refNet = exampleProject.getRunConfiguration(runConfigurationID).getReferenceNetwork();
		Collection<ConditionNetwork> condNets = exampleProject.getRunConfiguration(runConfigurationID).getConditionNetworks();
		
		CyNetworkManager networkManager = services.getCyNetworkManager();
		CyRootNetworkManager rootNetworkManager = services.getCyRootNetworkManager();
		CyNetworkViewFactory viewFactory = services.getCyNetworkViewFactory();
		CyNetworkViewManager viewManager = services.getCyNetworkViewManager();
		CyLayoutAlgorithm layout = services.getCyLayoutAlgorithmManager().getLayout("force-directed");
		TaskManager<?, ?> tm = services.getTaskManager();
		
		CyNetwork cyRefNet = CyNetworkBridge.createCyNetwork(refNet, services.getCyNetworkFactory());
		CyRootNetwork collection = rootNetworkManager.getRootNetwork(cyRefNet);
		collection.getRow(collection).set(CyRootNetwork.NAME, exampleProject.getName());
		networkManager.addNetwork(cyRefNet);
		
		CyNetworkView refView = viewFactory.createNetworkView(cyRefNet);
		viewManager.addNetworkView(refView);
		TaskIterator it = layout.createTaskIterator(refView, layout.createLayoutContext(), 
				CyLayoutAlgorithm.ALL_NODE_VIEWS, null);
		tm.execute(it);
		
		for (ConditionNetwork condNet : condNets){
			CyNetwork cyCondNet = CyNetworkBridge.createCyNetwork(condNet, collection);
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
