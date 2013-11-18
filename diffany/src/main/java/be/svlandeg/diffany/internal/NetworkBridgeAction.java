package be.svlandeg.diffany.internal;

import java.awt.event.ActionEvent;

import org.cytoscape.application.swing.AbstractCyAction;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.view.model.CyNetworkView;

import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.cytoscape.CyNetworkBridge;

public class NetworkBridgeAction extends AbstractCyAction {

	private Network network;
	private Services services;

	public NetworkBridgeAction(Services services, String menuTitle, Network network) {
		super(menuTitle, services.getCyApplicationManager(), null, null);
		setPreferredMenu("Apps.Diffany");
		this.services = services;
		this.network = network;
	}

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	@Override
	public void actionPerformed(ActionEvent e) {
		CyNetworkBridge bridge = new CyNetworkBridge(services.getCyNetworkFactory());
		CyNetwork cyNetwork = bridge.createCyNetwork(network);
		services.getCyNetworkManager().addNetwork(cyNetwork);
		CyNetworkView cyView = services.getCyNetworkViewFactory().createNetworkView(cyNetwork);
		services.getCyNetworkViewManager().addNetworkView(cyView);
	}

}
