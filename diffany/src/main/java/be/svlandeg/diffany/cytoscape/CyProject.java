package be.svlandeg.diffany.cytoscape;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CySubNetwork;
import org.cytoscape.view.layout.CyLayoutAlgorithm;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskManager;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.OverlappingNetwork;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.concepts.RunConfiguration;
import be.svlandeg.diffany.cytoscape.internal.Services;
import be.svlandeg.diffany.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.semantics.DefaultNodeMapper;



/**
 * TODO: update documentation
 * 
 * @author Thomas Van Parys
 *
 */
public class CyProject{
	
	/**
	 * Internal class to link a differential network to its overlap counterpart.
	 * 
	 * @author Thomas Van Parys
	 *
	 */
	public class CyNetworkPair{
		public CyNetworkPair(CyNetwork diffNet, CyNetwork overlapNet) {
			this.diffNet = diffNet;
			this.overlapNet = overlapNet;
		}
		public CyNetwork diffNet;
		public CyNetwork overlapNet;
	}
	
	private CyNetwork referenceNetwork;
	private Set<CyNetwork> conditionalNetworks = new HashSet<CyNetwork>();
	private Set<CyNetworkPair> resultNetworks = new HashSet<CyNetworkPair>();

	public static final String DEFAULT_PROJECT_NAME = "New Project";
	
	/**
	 * The Project used to run the algorithm and keep track of the ontologies
	 */
	Project project;
	/**
	 * A project keeps tracks of information about one collection of {@link CyNetwork}s
	 */
	private CyRootNetwork collection;
	private int latestRunConfigID;

	
	/**
	 * Construct empty project, named after the collection
	 */
	public CyProject(CyRootNetwork collection){
		this.collection = collection;
		String name = collection.getRow(collection).get(CyNetwork.NAME, String.class);
		this.project = new Project(name, new DefaultEdgeOntology(), new DefaultNodeMapper());
	}
	
	/**
	 * Gets the reference network for this project.
	 * 
	 * @return the reference network for this project. Returns null if no reference has been selected yet.
	 */
	public CyNetwork getReferenceNetwork() {
		return referenceNetwork;
	}

	/**
	 * Set the referenceNetwork for this project
	 * @param referenceNetwork
	 */
	public void setReferenceNetwork(CyNetwork referenceNetwork) {
		this.referenceNetwork = referenceNetwork;
	}

	/**
	 * Remove a conditional network from the list
	 * @param cyNetwork network to be removed
	 */
	public void removeConditionalNetwork(CyNetwork cyNetwork){
		this.conditionalNetworks.remove(cyNetwork);
	}


	/**
	 * Get the set of conditional networks for this project
	 * @return the set of conditional networks for this project
	 */
	public Set<CyNetwork> getConditionalNetworks() {
		return conditionalNetworks;
	}

	/**
	 * Set the set of conditional networks for this project
	 * @param conditionalNetworks set of conditional networks for this project
	 */
	public void setConditionalNetworks(Set<CyNetwork> conditionalNetworks) {
		this.conditionalNetworks = conditionalNetworks;
	}

	/**
	 * Add a conditional network to this {@link CyProject}
	 * @param conditionalNetwork
	 */
	public void addConditionalNetwork(CyNetwork conditionalNetwork){
		this.conditionalNetworks.add(conditionalNetwork);
	}
	
	private void addResultPair(CyNetwork cyDiffNet, CyNetwork cyOverlapNet) {
		this.resultNetworks.add(new CyNetworkPair(cyDiffNet, cyOverlapNet));
	}
	
	/**
	 * Generates a {@link RunConfiguration} and adds it to the {@link Project}
	 * 
	 * @return the ID of the generated {@link RunConfiguration}
	 * @throws InvalidRunConfigurationException is thrown when not all necessary parameters are there to construct the {@link RunConfiguration}
	 */
	public int generateRunConfiguration() throws InvalidRunConfigurationException{
		if (!canExecute()){
			throw new InvalidRunConfigurationException();
		}
		
		CyNetworkBridge bridge = new CyNetworkBridge();
		ReferenceNetwork refNet = bridge.getReferenceNetwork(this.getReferenceNetwork(), 
				project.getEdgeOntology(), project.getNodeMapper());
		
		Set<ConditionNetwork> condNets = new HashSet<ConditionNetwork>();
		for (CyNetwork cyCondNet : this.getConditionalNetworks()){
			ConditionNetwork condNet = bridge.getConditionNetwork(cyCondNet, 
					project.getEdgeOntology(), project.getNodeMapper());
			condNets.add(condNet);			
		}
		
		RunConfiguration rc = new RunConfiguration (refNet, condNets);
		int runConfigID = project.addRunConfiguration(rc);
		
		this.latestRunConfigID = runConfigID;
		
		return runConfigID;
	}

	

	

	/**
	 * Checks if all necessary parameters are set to execute the algorithm.
	 * 
	 * @return true when minimal input parameters are set
	 */
	public boolean canExecute(){
		return this.referenceNetwork!=null && !this.conditionalNetworks.isEmpty();
	}

	
	/**
	 * Get the set of resulting network pairs.
	 * @return the set of resulting network pairs.
	 */
	public Set<CyNetworkPair> getResultNetworks() {
		return resultNetworks;
	}
	
	/**
	 * Iterates over all {@link CyNetwork}s and creates a {@link Set} of all the used
	 * interactions (both from source as differential networks). This set can be used to define the visual mappings. 
	 * @return a Set of all used interactions in this project.
	 */
	public Set<String> getAllInteractions(){
		Set<String> interactions = new HashSet<String>();
		for (CyNetwork net : getAllNetworks()){
			interactions.addAll(this.getInteractions(net));
		}
		return interactions;
	}
	
	/**
	 * Get all interactions within a single network
	 *   
	 * @param net a {@link CyNetwork} contained in this {@link CyProject}
	 * @return a Set of all used interactions in the given {@link CyNetwork}
	 */
	private Collection<String> getInteractions(CyNetwork net) {
		Set<String> interactions = new HashSet<String>();
		for (CyEdge edge : net.getEdgeList()){
			String type = net.getRow(edge).get(CyEdge.INTERACTION, String.class);
			interactions.add(type);
		}
		return interactions;
	}


	/**
	 * Returns all {@link CyNetwork}s in this project. This includes the reference network, 
	 * the list of conditional networks and the generated (or loaded) differential networks and their
	 * overlap counterparts.
	 * 
	 * @return all {@link CyNetwork}s in this project, regardless of type.
	 */
	private Set<CyNetwork> getAllNetworks(){
		Set<CyNetwork> networks = new HashSet<CyNetwork>();
		if (this.referenceNetwork !=null){
			networks.add(this.referenceNetwork);			
		}
		for (CyNetwork conNet : this.conditionalNetworks){
			networks.add(conNet);
		}
		for (CyNetworkPair netPair : this.resultNetworks){
			networks.add(netPair.diffNet);
			networks.add(netPair.overlapNet);
		}
		return networks;
	}
	
	/**
	 * Returns all {@link CyNetworkView}s that correspond to the source networks (reference, conditional and overlap) in this project.
	 * @param services the collection of Cytoscape services.
	 * @return all {@link CyNetworkView}s that correspond to the networks in this project.
	 */
	public Set<CyNetworkView> getAllSourceViews(Services services){
		Set<CyNetworkView> views = new HashSet<CyNetworkView>();

		CyNetworkViewManager viewManager = services.getCyNetworkViewManager();
		if (this.referenceNetwork !=null){
			Collection<CyNetworkView> refViews = viewManager.getNetworkViews(referenceNetwork);
			views.addAll(refViews);			
		}
		
		for (CyNetwork condNet : conditionalNetworks){
			views.addAll(viewManager.getNetworkViews(condNet));
		}
		for (CyNetworkPair resPair : resultNetworks){
			views.addAll(viewManager.getNetworkViews(resPair.overlapNet));
		}
		return views;
	}
	
	/**
	 * Returns all {@link CyNetworkView}s that correspond to differential networks in this project.
	 * @param services the collection of Cytoscape services.
	 * @return all {@link CyNetworkView}s that correspond to the networks in this project.
	 */
	public Set<CyNetworkView> getAllDifferentialViews(Services services){
		Set<CyNetworkView> views = new HashSet<CyNetworkView>();

		CyNetworkViewManager viewManager = services.getCyNetworkViewManager();
		
		for (CyNetworkPair resPair : resultNetworks){
			views.addAll(viewManager.getNetworkViews(resPair.diffNet));
		}
		return views;
	}

	
	/**
	 * Checks the default edge and node tables of the networks in this project for a given SUID.
	 * 
	 * @return true if one of the tables in use has the given SUID
	 */
	public boolean containsTableId(long suid){
		Set<Long> idSet = new HashSet<Long>();
		for (CyNetwork net : this.getAllNetworks()){
			idSet.add(net.getDefaultEdgeTable().getSUID());
			idSet.add(net.getDefaultNodeTable().getSUID());
		}
		return idSet.contains(suid);
	}


	/**
	 * Iterates the {@link CyNetwork}s in the project and removes those that 
	 * have been destroyed in the Cytoscape session.
	 */
	public void removeDestroyedNetworks() {
		// TODO delete networks from the project when destroyed
		
	}

	/**
	 * 
	 * @return network collection containing the {@link CyNetwork}s for this project.
	 */
	public CyRootNetwork getCollection() {
		return collection;
	}

	/**
	 * @return the ID of the {@link RunConfiguration} that was generated and added to the {@link Project} 
	 * the last.
	 */
	public int getLatestRunConfigID() {
		return latestRunConfigID;
	}
	
	public String getName(){
		return this.project.getName();
	}

	@Override
	public String toString() {
		System.out.println("Loading name: " + this.getName());
		return this.getName();
	}

	/**
	 * After running a {@link RunConfiguration}, the newly generated differential and overlap networks
	 * should be transformed into {@link CyNetwork}s and added to this {@link CyProct}
	 */
	public void update(Services services) {
		RunConfiguration runConfig = project.getRunConfiguration(latestRunConfigID);
		this.addDifferentialNetworks(runConfig.getDifferentialNetworks(), services);
	}
	
	
	/**
	 * Convert and add the resulting networks after running the algorithm. Add these networks to the set of 
	 * results of this {@link CyProject}
	 * 
	 * @param differentialNetworks
	 * @param cyProject
	 */
	private void addDifferentialNetworks(Collection<DifferentialNetwork> differentialNetworks, Services services) {
		CyNetworkBridge bridge = new CyNetworkBridge();
		for (DifferentialNetwork network : differentialNetworks){
			
			//add the diffnet
			CyNetwork cyDiffNet = this.addCyNetwork(bridge, network, services);
			
			//add the overlap
			OverlappingNetwork overlap = network.getOverlappingNetwork();
			CyNetwork cyOverlapNet = this.addCyNetwork(bridge, overlap, services);
			
			this.addResultPair(cyDiffNet, cyOverlapNet);
		}
		
	}
	
	

	/**
	 * Converts a {@link Network} to a {@link CyNetwork} and registers it
	 * with Cytoscape as a {@link CySubNetwork} of the currently selected Network Collection.
	 * 
	 * 
	 * @param bridge the {@link Network} to {@link CyNetwork} convertor
	 * @param network the {@link Network} to be added
	 * @return the created and added {@link CyNetwork}
	 */
	private CyNetwork addCyNetwork(CyNetworkBridge bridge, Network network, Services services){
		CyNetwork cyNet = bridge.createCyNetwork(network, collection);
		
		services.getCyNetworkManager().addNetwork(cyNet);
		
		CyNetworkView cyView = services.getCyNetworkViewFactory().createNetworkView(cyNet);
		services.getCyNetworkViewManager().addNetworkView(cyView);
		
		CyLayoutAlgorithm layout = services.getCyLayoutAlgorithmManager().getLayout("force-directed");
		TaskIterator it = layout.createTaskIterator(cyView, layout.createLayoutContext(), CyLayoutAlgorithm.ALL_NODE_VIEWS, null);
		TaskManager<?, ?> tm = services.getTaskManager();
		tm.execute(it);
		
		return cyNet;
	}

	public Project getProject() {
		return project;
	}
	
	
	

	
	
}
