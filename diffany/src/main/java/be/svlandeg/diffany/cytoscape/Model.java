package be.svlandeg.diffany.cytoscape;

import java.util.Collection;
import java.util.HashSet;
import java.util.Observable;
import java.util.Set;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;
import org.cytoscape.view.model.CyNetworkView;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.cytoscape.gui.GUIModel;
import be.svlandeg.diffany.internal.CyActivator;
import be.svlandeg.diffany.internal.Services;
import be.svlandeg.diffany.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.semantics.DefaultNodeMapper;

/**
 * Model that keeps track of all settings and selections within the Cytoscape App.
 * Only calling {@link runAlgorithm} should construct the actual model to do the 
 * calculations and produce results, which are handed back to this model.
 * 
 * @author thpar
 *
 */
public class Model extends Observable{

	Services services;
	private GUIModel guiModel;
	
	private Project currentProject;
	
	public Model(Services services){
		this.services = services;
		this.guiModel = new GUIModel();
	}
	
	/**
	 * 
	 * @return the guiModel for this app.
	 */
	public GUIModel getGuiModel(){
		return this.guiModel;
	}
	
	/** 
	 * Use information from this model to construct the {@link Project} and 
	 * {@link Network}s to run the actual algorithm. 
	 */
	public void runAlgorithm(){
		CyNetworkBridge bridge = new CyNetworkBridge(this);
		CyNetwork cyRefNetwork = guiModel.getReferenceNetwork();
		ReferenceNetwork refNet = bridge.getReferenceNetwork(cyRefNetwork);
		
		System.out.println("Refnet");
		System.out.println(refNet.getStringRepresentation());
		System.out.println(refNet.writeEdgesTab());
		
		Set<CyNetwork> cyCondNetworks = guiModel.getConditionEntries();
		
		System.out.println("Condnets");
		Set<ConditionNetwork> condSet = new HashSet<ConditionNetwork>();
		for (CyNetwork cyNet : cyCondNetworks){
			ConditionNetwork condNet = bridge.getConditionNetwork(cyNet);
			condSet.add(condNet);			
			System.out.println(condNet.getStringRepresentation());
			System.out.println(condNet.writeEdgesTab());
			System.out.println("---");
		}
		
		
		
		currentProject = new Project("Default Project", refNet, condSet, new DefaultEdgeOntology(), new DefaultNodeMapper());
		
		double cutoff = 0.25;
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(currentProject, cutoff);
		
		addDifferentialNetworks(currentProject.getDifferentialNetworks());
	}
	
		
	private void addDifferentialNetworks(Collection<DifferentialNetwork> differentialNetworks) {
		System.out.println("Diffnets");
		CyNetworkBridge bridge = new CyNetworkBridge(this);
		for (Network network : differentialNetworks){
			System.out.println(network.getStringRepresentation());
			System.out.println(network.writeEdgesTab());
			CyNetwork cyNet = bridge.createCyNetwork(network);
			services.getCyNetworkManager().addNetwork(cyNet);
			CyNetworkView cyView = services.getCyNetworkViewFactory().createNetworkView(cyNet);
			services.getCyNetworkViewManager().addNetworkView(cyView);
		}
	}

	/**
	 * Iterates over all networks in the current cytoscape session and returns the set of network collections.
	 * 
	 * @return the set of network collections.
	 */
	public Set<CyRootNetwork> getNetworkCollections(){
		Set<CyNetwork> allNetworks = services.getCyNetworkManager().getNetworkSet();
		CyRootNetworkManager rootNetManager = services.getCyRootNetworkManager();
		
		Set<CyRootNetwork> set = new HashSet<CyRootNetwork>();
		
		for (CyNetwork net : allNetworks){
			set.add(rootNetManager.getRootNetwork(net));
		}
		
		return set;
	}
	
	/**
	 * Returns the {@link CyNetwork} with given name.
	 * Returns null if no such network exists.
	 * 
	 * @param id the network name
	 * @return 
	 */
	public CyNetwork getNetworkByName(String id){
		Set<CyNetwork> allNetworks = services.getCyNetworkManager().getNetworkSet();
		for (CyNetwork net : allNetworks){
			String name = net.getRow(net).get(CyNetwork.NAME, String.class);
			if (name.equals(id)){
				return net;
			}
		}
		return null;
	}
	
	/**
	 * Gives access to a subset of the services offered in this context, as loaded in the {@link CyActivator}
	 * 
	 * @return
	 */
	public Services getServices() {
		return services;
	}

	
	
	
}
