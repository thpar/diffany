package be.svlandeg.diffany.cytoscape.gui;

import java.awt.Component;
import java.util.Observable;
import java.util.Observer;
import java.util.Set;

import javax.swing.Icon;
import javax.swing.JPanel;

import org.cytoscape.application.swing.CytoPanelComponent;
import org.cytoscape.application.swing.CytoPanelName;
import org.cytoscape.model.CyNetwork;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.internal.Services;

public class TabPane extends JPanel implements CytoPanelComponent, Observer{

	private static final long serialVersionUID = 1L;
	private Model model;
	private Services services;

	public TabPane(Model model){
		this.model = model;
		model.addObserver(this);
		services = this.model.getServices();
		

//		Set<CyNetwork> allNetworks = services.getCyNetworkManager().getNetworkSet();
//		for(CyNetwork net : allNetworks){
//			String netName = net.getRow(net).get(CyNetwork.NAME, String.class);
//			System.out.println(netName);
//		}

	}
	
	@Override
	public Component getComponent() {
		return this;
	}

	@Override
	public CytoPanelName getCytoPanelName() {
		return CytoPanelName.WEST;
	}

	@Override
	public String getTitle() {
		return new String("Diffany");
	}

	@Override
	public Icon getIcon() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void update(Observable o, Object arg) {
		// TODO Auto-generated method stub
		
	}
	
	
}
