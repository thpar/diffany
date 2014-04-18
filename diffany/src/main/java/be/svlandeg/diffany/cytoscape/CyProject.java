package be.svlandeg.diffany.cytoscape;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.CyNetworkViewManager;

import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.OverlappingNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.DifferentialOutput;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunConfiguration;
import be.svlandeg.diffany.core.project.RunDiffConfiguration;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.cytoscape.internal.Services;



/**
 * The CyProject keeps track of all project specific settings. It is linked to a Collection ({@link CyRootNetwork} at the one
 * hand, and one {@link Project} at the other. 
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
	
	//currently selected networks. These are used to update the VisualStyles and
	//to construct the run configuration.
	private CyNetwork referenceNetwork;
	private Set<CyNetwork> conditionalNetworks = new HashSet<CyNetwork>();
	private Set<CyNetworkPair> resultNetworks = new HashSet<CyNetworkPair>();

	
	/**
	 * The Project used to run the algorithm and keep track of the ontologies
	 */
	private Project project;
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
	 * Set the referenceNetwork for this project. Register it with the {@link Project}
	 * 
	 * @param referenceNetwork
	 */
	public void setReferenceNetwork(CyNetwork referenceNetwork) {
		this.referenceNetwork = referenceNetwork;
		if (referenceNetwork != null){
			this.registerReferenceNetwork(referenceNetwork);			
		}
	}

	/**
	 * Remove a conditional network from the list.
	 * 
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
	 * Set the set of conditional networks for this project. Register them with the {@link Project}
	 * @param conditionalNetworks set of conditional networks for this project
	 */
	public void setConditionalNetworks(Set<CyNetwork> conditionalNetworks) {
		this.conditionalNetworks = conditionalNetworks;
		for (CyNetwork condNet : conditionalNetworks){
			this.registerConditionNetwork(condNet);
		}
	}

	/**
	 * Add a resulting differential network, paired with an overlap network to this project.
	 * Registration with the {@link Project} is not needed as this should already be aware of the structure
	 * of output networks.
	 * 
	 * @param cyDiffNet
	 * @param cyOverlapNet
	 */
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
		
		ReferenceNetwork refNet = CyNetworkBridge.getReferenceNetwork(this.getReferenceNetwork(), 
				project.getEdgeOntology(), project.getNodeMapper());
		
		Set<ConditionNetwork> condNets = new HashSet<ConditionNetwork>();
		for (CyNetwork cyCondNet : this.getConditionalNetworks()){
			ConditionNetwork condNet = CyNetworkBridge.getConditionNetwork(cyCondNet, 
					project.getEdgeOntology(), project.getNodeMapper());
			condNets.add(condNet);			
		}
		
		int runConfigID = project.addRunConfiguration(refNet, condNets);
		
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
	 * interactions, both from source as differential networks.
	 * @return a Set of all used interactions in this project.
	 * 
	 * @deprecated There is no use case where you need all interactions. 
	 * Use getAllSourceInteractions() or getAllDifferentialInteractions() instead.
	 */
	public Set<String> getAllInteractions(){
		Set<String> interactions = new HashSet<String>();
		for (CyNetwork net : getAllNetworks()){
			interactions.addAll(this.getInteractions(net));
		}
		return interactions;
	}
	
	/**
	 * Gathers the set of interactions (edge types) for all currently selected source networks.
	 * @return a set of interactions currently used in the source networks of this project
	 */
	public Set<String> getAllSourceInteractions(){
		Set<String> interactions = new HashSet<String>();
		for (CyNetwork net : getAllSourceNetworks()){
			interactions.addAll(this.getInteractions(net));
		}
		return interactions;
	}
	
	/**
	 * Gathers the set of interactions (edge types) for all currently selected differential networks.
	 * @return a set of interactions currently used in the differential networks of this project
	 */
	public Set<String> getAllDifferentialInteractions(){
		Set<String> interactions = new HashSet<String>();
		for (CyNetwork net : getAllDifferentialNetworks()){
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
		//TODO cache the interactions instead of iterating all edges every time
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
	 * Returns all {@link CyNetwork}s in this project currently used as source. 
	 * This includes the reference network, the list of conditional networks and the generated (or loaded) overlap networks
	 *  
	 * @return all {@link CyNetwork}s in this project, used as source.
	 */
	private Set<CyNetwork> getAllSourceNetworks(){
		Set<CyNetwork> networks = new HashSet<CyNetwork>();
		if (this.referenceNetwork !=null){
			networks.add(this.referenceNetwork);			
		}
		for (CyNetwork conNet : this.conditionalNetworks){
			networks.add(conNet);
		}
		for (CyNetworkPair netPair : this.resultNetworks){
			networks.add(netPair.overlapNet);
		}
		return networks;
	}
	
	/**
	 * Returns all {@link CyNetwork}s in this project currently used as differential network.
	 *  
	 * @return all differential {@link CyNetwork}s in this project.
	 */
	private Set<CyNetwork> getAllDifferentialNetworks(){
		Set<CyNetwork> networks = new HashSet<CyNetwork>();
		for (CyNetworkPair netPair : this.resultNetworks){
			networks.add(netPair.diffNet);
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
		
		for (CyNetwork net : this.getAllSourceNetworks()){
			views.addAll(viewManager.getNetworkViews(net));
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
		
		for (CyNetwork net : this.getAllDifferentialNetworks()){
			views.addAll(viewManager.getNetworkViews(net));
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
	 * Returns the network in this project with default (edge or node) table with given suid.
	 * 
	 * @param suid the id of a default edge or node table
	 * @return the network with the default edge or node table with given suid. 
	 * Returns null if the table is not used in this project.
	 */
	public CyNetwork getNetworkByTableId(Long suid) {
		CyNetwork foundNet = null;
		for (CyNetwork net : this.getAllNetworks()){
			if (net.getDefaultEdgeTable().getSUID() == suid || 
					net.getDefaultNodeTable().getSUID() == suid){
				foundNet = net;
			}
		}
		return foundNet;
	}
	
	/**
	 * @param net {@link CyNetwork} to be checked.
	 * @return true if the given {@link CyNetwork} is used as source network 
	 */
	public boolean isSourceNetwork(CyNetwork net){
		return this.getAllSourceNetworks().contains(net);
	}
	
	/**
	 * @param net {@link CyNetwork} to be checked.
	 * @return true if the given {@link CyNetwork} is used as differential network 
	 */

	public boolean isDiffNetwork(CyNetwork net){
		return this.getAllDifferentialNetworks().contains(net);
	}
	
	/**
	 * @param net {@link CyNetwork} to be checked.
	 * @return true if the given {@link CyNetwork} is the reference network 
	 */
	public boolean isReferenceNetwork(CyNetwork net){
		return net == this.referenceNetwork;
	}
	
	/**
	 * @param net {@link CyNetwork} to be checked.
	 * @return true if the given {@link CyNetwork} is a conditional network 
	 */
	public boolean isConditionalNetwork(CyNetwork net){
		return this.conditionalNetworks.contains(net);
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
	 * @return network collection (@link {@link CyRootNetwork}) containing the {@link CyNetwork}s for this project.
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
	
	/**
	 * Returns the name of this project. Delegated to the name of the {@link Project}.
	 * @return
	 */
	public String getName(){
		return this.project.getName();
	}

	@Override
	public String toString() {
		return this.getName();
	}

	/**
	 * Updates this {@link CyProject} after a run.
	 * After running a {@link RunConfiguration}, the newly generated differential and overlap networks
	 * should be transformed into {@link CyNetwork}s and added to this {@link CyProct}
	 */
	public void update(Services services) {
		// TODO (after Sofie's refactoring) check type of current ronConfigurationID
		RunDiffConfiguration runConfig = (RunDiffConfiguration) project.getRunConfiguration(latestRunConfigID);
		DifferentialOutput differentialOutput = runConfig.getDifferentialOutput();
		this.addDifferentialNetworks(differentialOutput, services);
	}
	
	
	/**
	 * Convert and add the resulting networks after running the algorithm. 
	 * Add these networks to the set of results of this {@link CyProject}
	 * 
	 * @param differentialOutput
	 * @param services 
	 */
	private void addDifferentialNetworks(DifferentialOutput differentialOutput, Services services) {
		for (OutputNetworkPair pair: differentialOutput.getOutputAsPairs()){
			DifferentialNetwork differentialNetwork = pair.getDifferentialNetwork();
			OverlappingNetwork overlappingNetwork = pair.getOverlappingNetwork();
			//add the diffnet
			CyNetwork cyDiffNet = CyNetworkBridge.addCyNetwork(differentialNetwork, this.collection, services);
				
			//add the overlap
			CyNetwork cyOverlapNet = CyNetworkBridge.addCyNetwork(overlappingNetwork, this.collection, services);
				
			this.addResultPair(cyDiffNet, cyOverlapNet);
		}
	}
	
	/**
	 * The actual {@link Project} that executes the algorithm.
	 * 
	 * @return The actual {@link Project} that executes the algorithm.
	 */
	public Project getProject() {
		return project;
	}
	
	/**
	 * When a reference network gets added or is altered in any way, the {@link Project} needs
	 * to be notified of this, in order to be able to update its ontologies.
	 * 
	 * @param network the {@link CyNetwork} containing new information
	 */
	public void registerReferenceNetwork(CyNetwork network){
		ReferenceNetwork net = CyNetworkBridge.getReferenceNetwork(network, this.project.getEdgeOntology(), this.project.getNodeMapper());
		
		// The logger object is needed as argument, but it's content will not be shown
		this.project.registerSourceNetwork(net, new Logger());
	}
	/**
	 * When a condition network gets added or is altered in any way, the {@link Project} needs
	 * to be notified of this, in order to be able to update its ontologies.
	 * 
	 * @param network the {@link CyNetwork} containing new information
	 */
	public void registerConditionNetwork(CyNetwork network){
		ConditionNetwork net = CyNetworkBridge.getConditionNetwork(network, this.project.getEdgeOntology(), this.project.getNodeMapper());
		
		// The logger object is needed as argument, but it's content will not be shown
		this.project.registerSourceNetwork(net, new Logger());
	}

	
	
	
}
