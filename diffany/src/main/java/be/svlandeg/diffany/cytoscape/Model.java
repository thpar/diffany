package be.svlandeg.diffany.cytoscape;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Observable;
import java.util.Set;

import javax.swing.JFrame;

import org.cytoscape.application.events.SetCurrentNetworkViewEvent;
import org.cytoscape.application.events.SetCurrentNetworkViewListener;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.events.NetworkAddedEvent;
import org.cytoscape.model.events.NetworkAddedListener;
import org.cytoscape.model.events.NetworkDestroyedEvent;
import org.cytoscape.model.events.NetworkDestroyedListener;
import org.cytoscape.model.events.RowsSetEvent;
import org.cytoscape.model.events.RowsSetListener;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.swing.DialogTaskManager;

import be.svlandeg.diffany.cytoscape.internal.CyActivator;
import be.svlandeg.diffany.cytoscape.internal.Services;
import be.svlandeg.diffany.cytoscape.tasks.UpdateVisualStyleTaskFactory;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualDiffStyle;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualSourceStyle;

/**
 * Model that keeps track of all settings and selections within the Cytoscape App.
 * It also acts as listener for Cytoscape events.
 *  
 * @author Thomas Van Parys
 *
 */
public class Model extends Observable implements NetworkAddedListener, 
												 NetworkDestroyedListener,
												 RowsSetListener,
												 SetCurrentNetworkViewListener{

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
	
	private boolean generateDiffNets = true;
	private boolean generateOverlapNets = true;

	private double cutoff = 0;

	private CyNetwork networkInFocus;


	
	/**
	 * Construct a new model and adds the app services to it. 
	 * Visual styles are initially registered with the VizMapper at this point.
	 * @param services the app services
	 */
	public Model(Services services){
		this.services = services;
	
		sourceStyle = new VisualSourceStyle(services);
		diffStyle = new VisualDiffStyle(services);
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
	 * Set the current {@link CyProject}. 
	 * 
	 * This change will trigger an update with all observers.
	 * 
	 * @param selectedProject 
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

	/**
	 * Get a list of all loaded {@link CyProject}s
	 * @return a list of all loaded {@link CyProject}s
	 */
	public Collection<CyProject> getProjects() {
		return projects.values();
	}
	
	/**
	 * Add a new {@link CyProject}. This happens only when a new network is added
	 * to the Cytoscape session as a new collection.
	 * 
	 * @param project
	 */
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

	
	/**
	 * Set the execution mode for the next run.
	 * @param mode {@link ComparisonMode}
	 */
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
	

	public boolean isGenerateDiffNets() {
		return generateDiffNets;
	}

	public void setGenerateDiffNets(boolean generateDiffNets) {
		this.generateDiffNets = generateDiffNets;
		setChanged();
		notifyObservers();
	}



	public boolean isGenerateOverlapNets() {
		return generateOverlapNets;
	}



	public void setGenerateOverlapNets(boolean generateOverlapNets) {
		this.generateOverlapNets = generateOverlapNets;
		setChanged();
		notifyObservers();
	}



	@Override
	public void handleEvent(SetCurrentNetworkViewEvent e) {
		//triggered when a CyView gets selected
		CyNetworkView view = e.getNetworkView();
		CyNetwork net = view.getModel();
		CyRootNetwork collection = services.getCyRootNetworkManager().getRootNetwork(net);
		this.networkInFocus = net;
		this.setSelectedProject(this.projects.get(collection));
	}


	/**
	 * The network of which the view is in focus.
	 * 
	 * @return The network of which the view is in focus. Null if there is no network to focus.
	 */
	public CyNetwork getNetworkInFocus() {
		return networkInFocus;
	}
	


}

