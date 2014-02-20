package be.svlandeg.diffany.cytoscape;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
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

import be.svlandeg.diffany.cytoscape.internal.CyActivator;
import be.svlandeg.diffany.cytoscape.internal.Services;
import be.svlandeg.diffany.cytoscape.tasks.UpdateVisualStyleTaskFactory;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualDiffStyle;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualSourceStyle;

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
	 * The currently selected project
	 */
	private CyProject selectedProject;

	/**
	 * All projects, one for each network collection
	 */
	private Map<CyRootNetwork, CyProject> projects = new HashMap<CyRootNetwork, CyProject>();
	
	
	private VisualSourceStyle sourceStyle;
	private VisualDiffStyle diffStyle;

	private JFrame swingApplication;
	
	/**
	 * Allows for switching between algorithm execution modes.
	 * 
	 * @author Thomas Van Parys
	 *
	 */
	public enum ComparisonMode{
		/**
		 * The reference network gets compared against each condition network separately
		 */
		REF_PAIRWISE,
		/**
		 * The reference network gets compared against all condition networks at once
		 */
		REF_TO_ALL;
		
		@Override
		public String toString(){
			switch(this){
			default:
			case REF_PAIRWISE:
				return "Pairwise";
			case REF_TO_ALL:
				return "One to all";
			}
		}
	}
	private ComparisonMode mode = ComparisonMode.REF_PAIRWISE;

	private double cutoff = 0;

	

	
	
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
	 * @return the current {@link CyProject}. null when no project has been selected yet or no collections are loaded at all.
	 */
	public CyProject getSelectedProject() {
		return selectedProject;
	}
	

	/**
	 * Set the collection of networks (aka the {@link CyRootNetwork}) that will be
	 * used for the algorithm.
	 * 
	 * This change will trigger an update with all observers.
	 * 
	 * @param selectedCollection
	 */
	public void setSelectedProject(CyProject selectedProject) {
		this.selectedProject = selectedProject;
		setChanged();
		notifyObservers();
	}


	@Override
	public void handleEvent(NetworkAddedEvent e) {
		//triggered on network added
		
		//add new project if needed
		CyRootNetworkManager rootNetManager = services.getCyRootNetworkManager();
		Set<CyRootNetwork> collections = projects.keySet();
		CyRootNetwork rootNet = rootNetManager.getRootNetwork(e.getNetwork());
		if (!collections.contains(rootNet)){
			this.addProject(new CyProject(rootNet));
		}
		
		setChanged();
		notifyObservers();
	}
	
	@Override
	public void handleEvent(NetworkDestroyedEvent e) {
		//triggered on network destroyed
		//TODO: destroyed networks should be removed from their projects
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
		//if yes, refresh the visual styles and re-apply them on the views
		Long suid = e.getSource().getSUID();
		
		if (this.selectedProject.containsTableId(suid)){
			CyNetwork net = this.selectedProject.getNetworkByTableId(suid);
			if (this.selectedProject.isReferenceNetwork(net)){
				this.selectedProject.registerReferenceNetwork(net);
			} else if (this.selectedProject.isConditionalNetwork(net)){
				this.selectedProject.registerConditionNetwork(net);
			}
			UpdateVisualStyleTaskFactory tf = new UpdateVisualStyleTaskFactory(this, this.selectedProject);
			TaskIterator it = tf.createTaskIterator();
			DialogTaskManager dtm = services.getDialogTaskManager();
			dtm.execute(it);		
		}
		
	}


	public Collection<CyProject> getProjects() {
		return projects.values();
	}
	
	
	public void addProject(CyProject project){
		this.projects.put(project.getCollection(), project);
	}


	/**
	 * Returns the execution mode for this project: Reference against all conditional networks, one at a time or all at once.
	 * @return execution (or comparison) mode of the project
	 */
	public ComparisonMode getMode() {
		return mode;
	}

	

	public void setMode(ComparisonMode mode) {
		this.mode = mode;
	}


	/**
	 * Get the edge cutoff.
	 * @return the edge cutoff.
	 */
	public double getCutoff() {
		return cutoff;
	}
	
	/**
	 * Set the edge cutoff
	 * @param cutoff the edge cutoff
	 */
	public void setCutoff(double cutoff) {
		this.cutoff  = cutoff;
	}


}
