package be.svlandeg.diffany.cytoscape.internal;

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

import java.util.HashMap;
import java.util.Map;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;
import org.cytoscape.view.layout.CyLayoutAlgorithmManager;
import org.cytoscape.view.model.CyNetworkViewFactory;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.view.vizmap.VisualMappingFunctionFactory;
import org.cytoscape.view.vizmap.VisualMappingManager;
import org.cytoscape.view.vizmap.VisualStyleFactory;
import org.cytoscape.work.TaskManager;
import org.cytoscape.work.swing.DialogTaskManager;
import org.osgi.framework.BundleContext;

/**
 * Class used to gather all services from the {@link BundleContext} as loaded in the {@link CyActivator}
 * 
 * @author Thomas Van Parys
 *
 */
public class Services {
	private CyApplicationManager cyApplicationManager;
	private CyNetworkFactory cyNetworkFactory;
	private CyNetworkManager cyNetworkManager;
	private CyNetworkViewFactory cyNetworkViewFactory;
	private CyNetworkViewManager cyNetworkViewManager;
	private TaskManager<?, ?> taskManager;
	private CyRootNetworkManager cyRootNetworkManager;
	private DialogTaskManager dialogTaskManager;
	private VisualMappingManager visualMappingManager;
	private VisualStyleFactory visualStyleFactory;
	private Map<String,VisualMappingFunctionFactory> visualMappingFunctionFactories =  
			new HashMap<String,VisualMappingFunctionFactory>();
	private CyLayoutAlgorithmManager cyLayoutAlgorithmManager;
		
	
	public void setCyApplicationManager(CyApplicationManager cyApplicationManager) {
		this.cyApplicationManager = cyApplicationManager;
	}


	public CyApplicationManager getCyApplicationManager() {
		return cyApplicationManager;
	}


	public CyNetworkFactory getCyNetworkFactory() {
		return cyNetworkFactory;
	}


	public void setCyNetworkFactory(CyNetworkFactory cyNetworkFactory) {
		this.cyNetworkFactory = cyNetworkFactory;
	}


	public CyNetworkManager getCyNetworkManager() {
		return cyNetworkManager;
	}


	public void setCyNetworkManager(CyNetworkManager cyNetworkManager) {
		this.cyNetworkManager = cyNetworkManager;
	}


	public CyNetworkViewFactory getCyNetworkViewFactory() {
		return cyNetworkViewFactory;
	}


	public void setCyNetworkViewFactory(CyNetworkViewFactory cyNetworkViewFactory) {
		this.cyNetworkViewFactory = cyNetworkViewFactory;
	}


	public CyNetworkViewManager getCyNetworkViewManager() {
		return cyNetworkViewManager;
	}


	public void setCyNetworkViewManager(CyNetworkViewManager cyNetworkViewManager) {
		this.cyNetworkViewManager = cyNetworkViewManager;
	}


	public void setTaskManager(TaskManager<?, ?> taskManager) {
		this.taskManager = taskManager;
		
	}


	public TaskManager<?, ?> getTaskManager() {
		return taskManager;
	}


	public void setCyRootNetworkManager(CyRootNetworkManager cyRootNetworkManager) {
		this.cyRootNetworkManager = cyRootNetworkManager;
	}


	public CyRootNetworkManager getCyRootNetworkManager() {
		return cyRootNetworkManager;
	}


	public void setDialogTaskManager(DialogTaskManager dialogTaskManager) {
		this.dialogTaskManager = dialogTaskManager;
		
	}


	public DialogTaskManager getDialogTaskManager() {
		return dialogTaskManager;
	}


	public VisualMappingManager getVisualMappingManager() {
		return visualMappingManager;
	}


	public void setVisualMappingManager(VisualMappingManager visualMappingManager) {
		this.visualMappingManager = visualMappingManager;
	}


	public VisualStyleFactory getVisualStyleFactory() {
		return visualStyleFactory;
	}


	public void setVisualStyleFactory(VisualStyleFactory visualStyleFactory) {
		this.visualStyleFactory = visualStyleFactory;
	}


	public VisualMappingFunctionFactory getVisualMappingFunctionFactory(String type) {
		return visualMappingFunctionFactories.get(type);
	}


	public void putVisualMappingFunctionFactory(String type,
			VisualMappingFunctionFactory visualMappingFunctionFactory) {
		this.visualMappingFunctionFactories.put(type, visualMappingFunctionFactory);
	}


	public void setCyLayoutAlgorithmManager(CyLayoutAlgorithmManager service) {
		this.cyLayoutAlgorithmManager = service;
	}


	public CyLayoutAlgorithmManager getCyLayoutAlgorithmManager() {
		return cyLayoutAlgorithmManager;
	}
	
	
}
