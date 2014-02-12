package be.svlandeg.diffany.cytoscape;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.CyNetworkViewManager;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.concepts.RunConfiguration;
import be.svlandeg.diffany.cytoscape.internal.Services;
import be.svlandeg.diffany.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.semantics.DefaultNodeMapper;



/**
 * The CyProject keeps the GUI settings for the current project and is able to generate a 
 * {@link Project} object to execute the algorithm.
 * 
 * @author Thomas Van Parys
 *
 */
public class CyProject extends Project{
	
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


	/**
	 * Default project name
	 */
	public static final String DEFAULT_PROJECT_NAME = "New Project";
	
	private double cutoff = 0;
	
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
	
	/**
	 * Construct empty project
	 */
	public CyProject(){
		super(DEFAULT_PROJECT_NAME, new DefaultEdgeOntology(), new DefaultNodeMapper());
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
				this.edgeOntology, this.nodeMapper);
		
		Set<ConditionNetwork> condNets = new HashSet<ConditionNetwork>();
		for (CyNetwork cyCondNet : this.getConditionalNetworks()){
			ConditionNetwork condNet = bridge.getConditionNetwork(cyCondNet, 
					this.edgeOntology, this.nodeMapper);
			condNets.add(condNet);			
		}
		
		RunConfiguration rc = new RunConfiguration (refNet, condNets);
		int runConfigID = this.addRunConfiguration(rc);
				
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
	 * Add a pair of resulting networks to the {@link CyProject}, probably after the algoritm has been run.
	 * A result pair consists of a differential network and its overlap counterpart.
	 * 
	 * @param diffNet the resulting differential network
	 * @param overlapNet its overlapping counterpart
	 */
	public void addResultPair(CyNetwork diffNet, CyNetwork overlapNet) {
		this.resultNetworks.add(new CyNetworkPair(diffNet, overlapNet));
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
	 * Returns the execution mode for this project: Reference against all conditional networks, one at a time or all at once.
	 * @return execution (or comparison) mode of the project
	 */
	public ComparisonMode getMode() {
		return mode;
	}

	/**
	 * Set the comparison mode for this project.
	 * 
	 * @param mode
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
		this.cutoff = cutoff;
	}

	/**
	 * Iterates the {@link CyNetwork}s in the project and removes those that 
	 * have been destroyed in the Cytoscape session.
	 */
	public void removeDestroyedNetworks() {
		// TODO Auto-generated method stub
		
		
	}

	
	
}
