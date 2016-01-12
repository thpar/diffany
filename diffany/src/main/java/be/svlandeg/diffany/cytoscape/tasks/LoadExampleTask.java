package be.svlandeg.diffany.cytoscape.tasks;

/*
 * #%L
 * Diffany
 * %%
 * Copyright (C) 2014 PSB/UGent - Sofie Van Landeghem and Thomas Van Parys
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU General Lesser Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/lgpl-3.0.html>.
 * #L%
 */

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

import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunDiffConfiguration;
import be.svlandeg.diffany.cytoscape.CyNetworkBridge;
import be.svlandeg.diffany.cytoscape.internal.Services;
import be.svlandeg.diffany.examples.GenericExample;

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
	private GenericExample example;
	
	private String exampleName = "";
	
	/**
	 * Construct a new task to load networks from a {@link Project}. The resulting networks and 
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
	
	/**
	 * Construct a new task to load networks from a {@link Project}. The resulting networks and 
	 * project settings will be ignored. Only the source networks are loaded and constructed. 
	 * 
	 * @param services the app {@link Services}
	 * @param example {@link Project} to be used as example input.
	 */
	public LoadExampleTask(Services services, GenericExample example) {
		this.services = services;
		this.example = example;
		this.exampleName = example.getName();
	}
	
	@Override
	public void run(TaskMonitor taskMonitor) throws Exception {
		taskMonitor.setTitle("Loading example: "+this.exampleName);			
		if (exampleProject == null){
			taskMonitor.setStatusMessage("Loading example networks");
			this.example.setTaskMonitor(taskMonitor);
			this.exampleProject = example.getDefaultProject();
			this.runConfigurationID = example.getDefaultRunConfigurationID(this.exampleProject);
		}
		
		taskMonitor.setStatusMessage("Constructing Cytoscape networks");
		ReferenceNetwork refNet = ((RunDiffConfiguration) exampleProject.getRunConfiguration(runConfigurationID)).getReferenceNetwork();
		Collection<ConditionNetwork> condNets = ((RunDiffConfiguration) exampleProject.getRunConfiguration(runConfigurationID)).getConditionNetworks();
		
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
		//No cancel feature		
	}

}
