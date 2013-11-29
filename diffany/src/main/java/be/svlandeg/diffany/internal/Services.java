package be.svlandeg.diffany.internal;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;
import org.cytoscape.view.model.CyNetworkViewFactory;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.work.TaskManager;
import org.cytoscape.work.swing.DialogTaskManager;
import org.osgi.framework.BundleContext;

/**
 * Class used to gather all services from the {@link BundleContext} as loaded in the {@link CyActivator}
 * 
 * @author thpar
 *
 */
public class Services {
	private CyApplicationManager cyApplicationManager;
	private CyNetworkFactory cyNetworkFactory;
	private CyNetworkManager cyNetworkManager;
	private CyNetworkViewFactory cyNetworkViewFactory;
	private CyNetworkViewManager cyNetworkViewManager;
	private TaskManager taskManager;
	private CyRootNetworkManager cyRootNetworkManager;
	private DialogTaskManager dialogTaskManager;
	
	
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


	public void setTaskManager(TaskManager taskManager) {
		this.taskManager = taskManager;
		
	}


	public TaskManager getTaskManager() {
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
	
	
	
	
	
	
}
