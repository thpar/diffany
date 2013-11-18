package be.svlandeg.diffany.internal;

import java.util.Properties;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.service.util.AbstractCyActivator;
import org.cytoscape.view.model.CyNetworkViewFactory;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.work.TaskManager;
import org.osgi.framework.BundleContext;

import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Node;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.cytoscape.actions.MenuAction;
import be.svlandeg.diffany.cytoscape.actions.NetworkBridgeAction;

/**
 * Defines a MenuAction and registers it as an OSGi service to the Cytoscape
 * application manager.
 * 
 * Generated from the cyaction-app Archetype.
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
		services.setCyNetworkViewManager(getService(context, CyNetworkViewManager.class));
		services.setTaskManager(getService(context, TaskManager.class));
		
		//create and register testing menu
		MenuAction action = new MenuAction(services.getCyApplicationManager(), "Diffany");
		Properties properties = new Properties();
		registerAllServices(context, action, properties);
		
		//create a testing network
		//------
		Network testNetwork = new ReferenceNetwork("My testing network");
		Node a = new Node("A");
		Node b = new Node("B");
		Node c = new Node("C");
		Node d = new Node("D");
		Edge ab = new Edge("binds", a,b,true);
		Edge bc = new Edge("binds", b,c,true);
		Edge ad = new Edge("binds", a,d,true);
		testNetwork.addNode(a);
		testNetwork.addNode(b);
		testNetwork.addNode(c);
		testNetwork.addNode(d);
		testNetwork.addEdge(ab);
		testNetwork.addEdge(bc);
		testNetwork.addEdge(ad);
		//------
		
		//convert network to CyNetwork and register it and add it to the menu as well
		NetworkBridgeAction networkAction = new NetworkBridgeAction(services,"Network test", testNetwork);
		registerAllServices(context, networkAction, new Properties());
	}

}
