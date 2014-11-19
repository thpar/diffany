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
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.progress.ProgressListener;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunConfiguration;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
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
	 * Internal class to link a differential network to its consensus counterpart.
	 * A pair is not necessary complete. If only one of both has been generated, the counterpart will be null.
	 * 
	 * @author Thomas Van Parys
	 *
	 */
	public class CyNetworkPair{
		public CyNetworkPair(CyNetwork diffNet, CyNetwork consensusNet) {
			this.diffNet = diffNet;
			this.consensusNet = consensusNet;
		}
		public CyNetwork diffNet;
		public CyNetwork consensusNet;
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
		this.project = new Project(name, new DefaultEdgeOntology());
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
	 * Set the referenceNetwork for this project. Register it with the {@link Project}. 
	 * 
	 * @param referenceNetwork the reference network for this {@link CyProject}. NULL if no reference network is selected.
	 */
	public void setReferenceNetwork(CyNetwork referenceNetwork) {
		if (this.referenceNetwork != referenceNetwork){
			this.referenceNetwork = referenceNetwork;
			if (referenceNetwork != null){
				this.registerReferenceNetwork(referenceNetwork);			
			}			
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
		//which conditional networks are actually new and need to be registered?
		HashSet<CyNetwork> newConditionalNetworks = new HashSet<CyNetwork>();
		for (CyNetwork condNet : conditionalNetworks){
			if (!this.conditionalNetworks.contains(condNet)){
				newConditionalNetworks.add(condNet);
			}
		}
		
		//set the new set
		this.conditionalNetworks = conditionalNetworks;
		
		//register only the new ones
		for (CyNetwork condNet : newConditionalNetworks){
			this.registerConditionNetwork(condNet);
		}
	}

	/**
	 * Add a resulting differential network, paired with an consensus network to this project.
	 * Registration with the {@link Project} is not needed as this should already be aware of the structure
	 * of output networks.
	 * 
	 * @param cyDiffNet
	 * @param cyConsensusNet
	 */
	private void addResultPair(CyNetwork cyDiffNet, CyNetwork cyConsensusNet) {
		this.resultNetworks.add(new CyNetworkPair(cyDiffNet, cyConsensusNet));
	}
	
	/**
	 * Generates a {@link RunConfiguration} and adds it to the {@link Project}
	 * 
	 * @return the ID of the generated {@link RunConfiguration}
	 * @throws InvalidRunConfigurationException is thrown when not all necessary parameters are there to construct the {@link RunConfiguration}
	 */
	public int generateRunConfiguration(Model model, ProgressListener listener) throws InvalidRunConfigurationException{
		if (!canExecute(model)){
			throw new InvalidRunConfigurationException();
		}
		
		
		int runConfigID;
		if (model.isGenerateDiffNets()){
			//"classic" ref/cond configuration
			ReferenceNetwork refNet = CyNetworkBridge.getReferenceNetwork(this.getReferenceNetwork(), 
					project.getEdgeOntology());						
			Set<ConditionNetwork> condNets = new HashSet<ConditionNetwork>();
			for (CyNetwork cyCondNet : this.getConditionalNetworks()){
				ConditionNetwork condNet = CyNetworkBridge.getConditionNetwork(cyCondNet, 
						project.getEdgeOntology());
				condNets.add(condNet);			
			}
			runConfigID = project.addRunConfiguration(refNet, condNets, model.getOverlapSupportCutoff(), true, listener);
			

					
			
		} else if (model.isGenerateConsensusNets()) {
			Set<InputNetwork> inputNetworks = new HashSet<InputNetwork>();
			
			//model might still mark one network as reference, if one was selected already. Treat it as input network in this configuration
			if (this.getReferenceNetwork() != null){
				inputNetworks.add(CyNetworkBridge.getReferenceNetwork(this.getReferenceNetwork(), project.getEdgeOntology()));
			}
			//then add all the conditional networks
			for (CyNetwork cyCondNet : this.getConditionalNetworks()){
				InputNetwork condNet = CyNetworkBridge.getInputNetwork(cyCondNet, 
						project.getEdgeOntology());
				inputNetworks.add(condNet);			
			}
			runConfigID = project.addRunConfiguration(inputNetworks, model.getOverlapSupportCutoff(), model.isRefIncludedInOverlapSupportCutoff(), listener);
			
		} else {
			//if we're not generating diff or consensus networks, then what are we doing here?
			throw new InvalidRunConfigurationException();
		}

		
		this.latestRunConfigID = runConfigID;
		
		return runConfigID;
	}

	

	

	/**
	 * Checks if all necessary parameters are set to execute the algorithm.
	 * 
	 * @return true when minimal input parameters are set
	 */
	public boolean canExecute(Model model){
		if (model.isGenerateDiffNets()){
			return this.referenceNetwork!=null && !this.conditionalNetworks.isEmpty();			
		} else if (model.isGenerateConsensusNets()){
			int numberOfNetworks = this.conditionalNetworks.size();
			if (this.referenceNetwork!=null){
				numberOfNetworks+=1;
			}
			return numberOfNetworks>=2;
		} else {
			return false;
		}
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
	 * consensus counterparts.
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
			if (netPair.diffNet != null){
				networks.add(netPair.diffNet);				
			}
			if (netPair.consensusNet!=null){
				networks.add(netPair.consensusNet);				
			}
		}
		return networks;
	}
	
	/**
	 * Returns all {@link CyNetwork}s in this project currently used as source. 
	 * This includes the reference network, the list of conditional networks and the generated (or loaded) consensus networks
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
			if (netPair.consensusNet !=null){
				networks.add(netPair.consensusNet);				
			}
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
			if (netPair.diffNet != null){
				networks.add(netPair.diffNet);				
			}
		}
		return networks;
	}
	/**
	 * Returns all {@link CyNetworkView}s that correspond to the source networks (reference, conditional and consensus) in this project.
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
	 * After running a {@link RunConfiguration}, the newly generated differential and consensus networks
	 * should be transformed into {@link CyNetwork}s and added to this {@link CyProject}
	 */
	public void update(Services services) {
		RunOutput runOutput = project.getOutput(latestRunConfigID);
		this.addDifferentialNetworks(runOutput, services);
	}
	
	
	/**
	 * Convert and add the resulting networks after running the algorithm. 
	 * Add these networks to the set of results of this {@link CyProject}
	 * 
	 * @param runOutput
	 * @param services 
	 */
	private void addDifferentialNetworks(RunOutput runOutput, Services services) {
		Set<DifferentialNetwork> addedDiffNets = new HashSet<DifferentialNetwork>();
		Set<ConsensusNetwork> addedConsensusNets = new HashSet<ConsensusNetwork>();
		
		for (OutputNetworkPair pair: runOutput.getOutputAsPairs()){
			DifferentialNetwork differentialNetwork = pair.getDifferentialNetwork();
			ConsensusNetwork consensusNetwork = pair.getConsensusNetwork();
			addedDiffNets.add(differentialNetwork);
			addedConsensusNets.add(consensusNetwork);
			
			//add the diffnet
			CyNetwork cyDiffNet = null; 
			if (differentialNetwork != null){
				cyDiffNet = CyNetworkBridge.addCyNetwork(differentialNetwork, this.collection, services);				
			}
				
			//add the consensus net
			CyNetwork cyConsensusNet = null;
			if (consensusNetwork != null){
				cyConsensusNet = CyNetworkBridge.addCyNetwork(consensusNetwork, this.collection, services);				
			}
				
			this.addResultPair(cyDiffNet, cyConsensusNet);
		}
		
		//add differential networks that were not added as a pair
		for (DifferentialNetwork diffNet : runOutput.getDifferentialNetworks()){
			if (!addedDiffNets.contains(diffNet)){
				CyNetwork cyDiffNet = CyNetworkBridge.addCyNetwork(diffNet, this.collection, services);
				this.addResultPair(cyDiffNet, null);
			}
		}
		
		//add consensus networks that were not added as a pair
		for (ConsensusNetwork consensusNet : runOutput.getConsensusNetworks()){
			if (!addedConsensusNets.contains(consensusNet)){
				CyNetwork cyConsensusNet = CyNetworkBridge.addCyNetwork(consensusNet, this.collection, services);
				this.addResultPair(null, cyConsensusNet);
			}
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
		ReferenceNetwork net = CyNetworkBridge.getReferenceNetwork(network, this.project.getEdgeOntology());
		
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
		ConditionNetwork net = CyNetworkBridge.getConditionNetwork(network, this.project.getEdgeOntology());
		
		// The logger object is needed as argument, but it's content will not be shown
		this.project.registerSourceNetwork(net, new Logger());
	}

	/**
	 * 
	 * @return the number of selected condition and reference networks 
	 */
	public int getNumberOfInputNetworks(){
		int condNumber = this.conditionalNetworks.size();
		if (this.referenceNetwork != null){
			condNumber++;
		}
		return condNumber;
		
	}
	
	
}
