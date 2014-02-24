package be.svlandeg.diffany.cytoscape.internal;

import java.util.Properties;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.application.events.SetCurrentNetworkViewListener;
import org.cytoscape.application.swing.CySwingApplication;
import org.cytoscape.application.swing.CytoPanelComponent;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.model.events.NetworkAddedListener;
import org.cytoscape.model.events.NetworkDestroyedListener;
import org.cytoscape.model.events.RowsSetListener;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;
import org.cytoscape.service.util.AbstractCyActivator;
import org.cytoscape.session.events.SessionAboutToBeSavedListener;
import org.cytoscape.session.events.SessionLoadedListener;
import org.cytoscape.view.layout.CyLayoutAlgorithm;
import org.cytoscape.view.layout.CyLayoutAlgorithmManager;
import org.cytoscape.view.model.CyNetworkViewFactory;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.view.vizmap.VisualMappingFunctionFactory;
import org.cytoscape.view.vizmap.VisualMappingManager;
import org.cytoscape.view.vizmap.VisualStyleFactory;
import org.cytoscape.work.TaskManager;
import org.cytoscape.work.swing.DialogTaskManager;
import org.osgi.framework.BundleContext;

import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.SessionListener;
import be.svlandeg.diffany.cytoscape.actions.LoadExampleAction;
import be.svlandeg.diffany.cytoscape.actions.RunProjectAction;
import be.svlandeg.diffany.cytoscape.gui.TabPane;
import be.svlandeg.diffany.cytoscape.layout.CopyLayout;
import be.svlandeg.diffany.examples.Bandyopadhyay2010;
import be.svlandeg.diffany.examples.Ideker2011;

/**
 * Entry point for the Diffany Cytoscape 3 App. Here the necessary services are called and bundled into the
 * {@link Model} to be used throughout the app. The services this app offers are created here and registered
 * within the {@link BundleContext}.
 * 
 * @author Thomas Van Parys
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
		services.setVisualStyleFactory(getService(context, VisualStyleFactory.class));
		services.setVisualMappingManager(getService(context, VisualMappingManager.class));
		services.putVisualMappingFunctionFactory("continuous", 
				getService(context, VisualMappingFunctionFactory.class,"(mapping.type=continuous)"));
		services.putVisualMappingFunctionFactory("discrete", 
				getService(context, VisualMappingFunctionFactory.class,"(mapping.type=discrete)"));
		services.putVisualMappingFunctionFactory("passthrough", 
				getService(context, VisualMappingFunctionFactory.class,"(mapping.type=passthrough)"));

		services.setCyLayoutAlgorithmManager(getService(context, CyLayoutAlgorithmManager.class));
		
		CySwingApplication swingApplication = getService(context, CySwingApplication.class);
		
		Model model = new Model(services);
		
		model.setParentWindow(swingApplication.getJFrame());
		
		//Create and register the control panel
		TabPane sidePane = new TabPane(model);
		
		registerService(context,sidePane,CytoPanelComponent.class, new Properties());
		
		//register action to run the current Diffany project
		RunProjectAction runProjectAction = new RunProjectAction(model,"Run Diffany project");
		registerAllServices(context, runProjectAction, new Properties());
		
		//register actions to import the projects defined in the examples package
		Bandyopadhyay2010 example1 = new Bandyopadhyay2010();
		Project exampleProject1 = example1.getProjectFigure1C();
		registerAllServices(context, new LoadExampleAction(services,"Bandyopadhyay2010", 
				exampleProject1, example1.getTestConfiguration1C(exampleProject1)), 
				new Properties());
		
		Ideker2011 example2 = new Ideker2011();
		Project exampleProject2 = example2.getProjectFigure3A();
		registerAllServices(context, new LoadExampleAction(services,"Ideker2011", 
				exampleProject2, example2.getTestConfiguration3A(exampleProject2)), 
				new Properties());
		
		//Register network listeners
		registerService(context,model, NetworkAddedListener.class, new Properties());
		registerService(context,model, NetworkDestroyedListener.class, new Properties());
		registerService(context,model, RowsSetListener.class, new Properties());
		
		registerService(context,model, SetCurrentNetworkViewListener.class, new Properties());
		
		//Register session handlers
		registerService(context, new SessionListener(model),SessionLoadedListener.class, new Properties() );
		registerService(context, new SessionListener(model),SessionAboutToBeSavedListener.class, new Properties() );
		
		//Register layouting algorithm
		CopyLayout copyLayout = new CopyLayout(model);
		Properties myLayoutProps = new Properties();
		registerService(context,copyLayout,CyLayoutAlgorithm.class, myLayoutProps);
		
	}

}
