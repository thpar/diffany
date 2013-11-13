package be.svlandeg.diffany.internal;

import java.util.Properties;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.service.util.AbstractCyActivator;
import org.cytoscape.view.model.CyNetworkViewFactory;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.osgi.framework.BundleContext;

import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Node;
import be.svlandeg.diffany.concepts.ReferenceNetwork;

/**
 * Defines a MenuAction and registers it as an OSGi service to the Cytoscape
 * application manager.
 * 
 * Generated from the cyaction-app Archetype.
 */
public class CyActivator extends AbstractCyActivator
{

	@Override
	public void start(BundleContext context) throws Exception
	{

		CyApplicationManager cyApplicationManager = getService(context, CyApplicationManager.class);
		CyNetworkManager cyNetworkManager = getService(context, CyNetworkManager.class);
		CyNetworkViewManager cyNetworkViewManager = getService(context, CyNetworkViewManager.class);
		CyNetworkFactory cyNetworkFactory = getService(context, CyNetworkFactory.class);
		CyNetworkViewFactory cyNetworkViewFactory = getService(context, CyNetworkViewFactory.class);
		
		//create and register testing menu
		MenuAction action = new MenuAction(cyApplicationManager, "Diffany");
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
		
		//convert network to CyNetwork and register it
		NetworkBridgeAction networkAction = new NetworkBridgeAction(cyApplicationManager,
				cyNetworkFactory, cyNetworkManager, cyNetworkViewManager, cyNetworkViewFactory,"Network test", testNetwork);
		
		registerAllServices(context, networkAction, new Properties());
	}

}
