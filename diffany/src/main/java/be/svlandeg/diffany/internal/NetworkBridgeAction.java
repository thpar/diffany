package be.svlandeg.diffany.internal;

import java.awt.event.ActionEvent;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.application.swing.AbstractCyAction;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.CyNetworkViewFactory;
import org.cytoscape.view.model.CyNetworkViewManager;

import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.cytoscape.bridge.CyNetworkBridge;

public class NetworkBridgeAction extends AbstractCyAction {

	private CyNetworkFactory cyNetworkFactory;
	private Network network;
	private CyNetworkManager cyNetworkManager;
	private CyNetworkViewFactory cyNetworkViewFactory;
	private CyNetworkViewManager cyNetworkViewManager;

	public NetworkBridgeAction(CyApplicationManager applicationManager, CyNetworkFactory cyNetworkFactory, 
			CyNetworkManager cyNetworkManager, CyNetworkViewManager cyNetworkViewManager, 
			CyNetworkViewFactory cyNetworkViewFactory, 
			String menuTitle, Network network) {
		super(menuTitle, applicationManager, null, null);
		setPreferredMenu("Apps.Diffany");
		this.cyNetworkFactory = cyNetworkFactory;
		this.cyNetworkViewFactory = cyNetworkViewFactory;
		this.cyNetworkManager = cyNetworkManager;
		this.cyNetworkViewManager = cyNetworkViewManager;
		
		this.network = network;
	}

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	@Override
	public void actionPerformed(ActionEvent e) {
		CyNetwork cyNetwork = cyNetworkFactory.createNetwork();
		CyNetworkBridge bridge = new CyNetworkBridge();
		bridge.convertToCyNetwork(cyNetwork, network);
		cyNetworkManager.addNetwork(cyNetwork);
		CyNetworkView cyView = cyNetworkViewFactory.createNetworkView(cyNetwork);
		cyNetworkViewManager.addNetworkView(cyView);
	}

}
