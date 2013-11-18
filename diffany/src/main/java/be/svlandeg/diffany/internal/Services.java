package be.svlandeg.diffany.internal;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.view.model.CyNetworkViewFactory;
import org.cytoscape.view.model.CyNetworkViewManager;

public class Services {
	private CyApplicationManager cyApplicationManager;
	private CyNetworkFactory cyNetworkFactory;
	private CyNetworkManager cyNetworkManager;
	private CyNetworkViewFactory cyNetworkViewFactory;
	private CyNetworkViewManager cyNetworkViewManager;
	
	
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
	
	
	
	
}
