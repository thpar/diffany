package be.svlandeg.diffany.cytoscape.gui;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Observable;
import java.util.Observer;
import java.util.Set;

import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JPanel;

import org.cytoscape.application.swing.CytoPanelComponent;
import org.cytoscape.application.swing.CytoPanelName;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CyRootNetwork;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.internal.Services;

public class TabPane extends JPanel implements CytoPanelComponent, Observer, ActionListener{

	private static final long serialVersionUID = 1L;
	private Model model;
	private Services services;

	public TabPane(Model model){
		this.model = model;
		model.addObserver(this);
		services = this.model.getServices();
		this.setSize(new Dimension(200,100));
		
		JButton button = new JButton("Push me");
		button.addActionListener(this);
		this.add(button);
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

	@Override
	public void actionPerformed(ActionEvent e) {
		Set<CyNetwork> allNetworks = services.getCyNetworkManager().getNetworkSet();
		for(CyNetwork net : allNetworks){
			int subNets = 0;
			String netName = net.getRow(net).get(CyNetwork.NAME, String.class);
			if (net instanceof CyRootNetwork){
				CyRootNetwork rootNet = (CyRootNetwork)net;
				subNets = rootNet.getSubNetworkList().size();
			}
			System.out.println(netName + " " + subNets);
			
		}
		
	}
	
	
}
