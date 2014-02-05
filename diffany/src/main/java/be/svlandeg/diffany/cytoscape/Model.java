package be.svlandeg.diffany.cytoscape;

import java.util.HashSet;
import java.util.Observable;
import java.util.Set;

import javax.swing.JFrame;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.events.NetworkAddedEvent;
import org.cytoscape.model.events.NetworkAddedListener;
import org.cytoscape.model.events.NetworkDestroyedEvent;
import org.cytoscape.model.events.NetworkDestroyedListener;
import org.cytoscape.model.events.RowsSetEvent;
import org.cytoscape.model.events.RowsSetListener;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.swing.DialogTaskManager;

import be.svlandeg.diffany.cytoscape.tasks.UpdateVisualStyleTaskFactory;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualDiffStyle;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualSourceStyle;
import be.svlandeg.diffany.internal.CyActivator;
import be.svlandeg.diffany.internal.Services;

/**
 * Model that keeps track of all settings and selections within the Cytoscape App.
 *  
 * @author Thomas Van Parys
 *
 */
public class Model extends Observable implements NetworkAddedListener, 
												 NetworkDestroyedListener,
												 RowsSetListener{

	/**
	 * A collection of all Cytoscape services that were registered in the {@link CyActivator}
	 */
	private Services services;
	
	/**
	 * The current project
	 */
	private CyProject currentProject = new CyProject();

	/**
	 * Network collection that's currently selected
	 */
	private CyRootNetwork selectedCollection;
	
	private VisualSourceStyle sourceStyle;
	private VisualDiffStyle diffStyle;

	private JFrame swingApplication;
	
	/**
	 * Construct a new model and adds the app services to it
	 * @param services the app services
	 */
	public Model(Services services){
		this.services = services;
	
		sourceStyle = new VisualSourceStyle(services);
		diffStyle = new VisualDiffStyle(services);
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
	 * @return {@link CyNetwork} with given id
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

	
	/**
	 * Returns the current {@link CyProject}.
	 * @return the current {@link CyProject}
	 */
	public CyProject getCurrentProject() {
		return currentProject;
	}
	
	
	/**
	 * The collection of networks (aka the {@link CyRootNetwork}) that will be
	 * used for the algorithm. 
	 * @return The collection of networks (aka the {@link CyRootNetwork}) that will be used for the algorithm. Returns null
	 * if no collection has been selected (yet).
	 */
	public CyRootNetwork getSelectedCollection() {
		return selectedCollection;
	}

	/**
	 * Set the collection of networks (aka the {@link CyRootNetwork}) that will be
	 * used for the algorithm.
	 * 
	 * This change will trigger an update with all observers.
	 * 
	 * @param selectedCollection
	 */
	public void setSelectedCollection(CyRootNetwork selectedCollection) {
		this.selectedCollection = selectedCollection;
		setChanged();
		notifyObservers();
	}


	@Override
	public void handleEvent(NetworkAddedEvent e) {
		//triggered on network added
		setChanged();
		notifyObservers();
	}
	
	@Override
	public void handleEvent(NetworkDestroyedEvent e) {
		//triggered on network destroyed
		setChanged();
		notifyObservers();
	}
	
	/**
	 * Get the style applied to source networks
	 * @return the style applied to source networks
	 */
	public VisualSourceStyle getSourceStyle() {
		return sourceStyle;
	}

	/**
	 * Get the style applied to differential networks
	 * @return the style applied to differential networks
	 */
	public VisualDiffStyle getDiffStyle() {
		return diffStyle;
	}

	/**
	 * Set a reference to the Cytoscape main window
	 * @param jFrame a reference to the Cytoscape main window
	 */
	public void setParentWindow(JFrame jFrame) {
		this.swingApplication = jFrame;
	}
	
	/**
	 * Get a reference to the Cytoscape main window
	 * @return a reference to the Cytoscape main window
	 */
	public JFrame getParentWindow(){
		return this.swingApplication;
	}


	@Override
	public void handleEvent(RowsSetEvent e) {
		//triggered when one or more rows change in a CyTable
		
		//check if the row is part of a table we care about				
		//if yes, refresh the visual styles and reapply them on the views
		Long suid = e.getSource().getSUID();
		if (this.currentProject.containsTableId(suid)){
			UpdateVisualStyleTaskFactory tf = new UpdateVisualStyleTaskFactory(this, this.getCurrentProject());
			TaskIterator it = tf.createTaskIterator();
			DialogTaskManager dtm = services.getDialogTaskManager();
			dtm.execute(it);		
		}
		
	}

	
}
