package be.svlandeg.diffany.internal;

import java.util.Properties;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.application.swing.CytoPanelComponent;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.model.events.NetworkAddedListener;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;
import org.cytoscape.service.util.AbstractCyActivator;
import org.cytoscape.view.model.CyNetworkViewFactory;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.work.TaskManager;
import org.cytoscape.work.swing.DialogTaskManager;
import org.osgi.framework.BundleContext;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.actions.RunProjectAction;
import be.svlandeg.diffany.cytoscape.gui.TabPane;

/**
 * Entry point for the Diffany Cytoscape App. Here the necessary services are called and bundled into the
 * {@link Model} to be used throughout the app. The services this app offers are created here and registered
 * within the {@link BundleContext}.
 */
public class CyActivator extends AbstractCyActivator
{

	@Override
	public void start(BundleContext context) throws Exception{

		//load all needed services
		Services services = new Services();
		services.setCyApplicationManager(getService(context, CyApplicationManager.class));
		services.setCyNetworkFactory(getService(context, CyNetworkFactory.class));
		services.setCyNetworkViewFactory(getService(context, CyNetworkViewFactory.class));
		services.setCyNetworkManager(getService(context, CyNetworkManager.class));
		services.setCyRootNetworkManager(getService(context, CyRootNetworkManager.class));
		services.setCyNetworkViewManager(getService(context, CyNetworkViewManager.class));
		services.setTaskManager(getService(context, TaskManager.class));
		services.setDialogTaskManager(getService(context, DialogTaskManager.class));
		
		
		Model model = new Model(services);
		
		
		//register action to run  the current Diffany project
		RunProjectAction runProjectAction = new RunProjectAction(model,"Run Diffany project");
		registerAllServices(context, runProjectAction, new Properties());
		
		
		//Create and register the control panel
		TabPane sidePane = new TabPane(model);
		registerService(context,sidePane,CytoPanelComponent.class, new Properties());
		
		//Register control panel as network listener
		registerService(context,sidePane, NetworkAddedListener.class, new Properties());
	}

}
