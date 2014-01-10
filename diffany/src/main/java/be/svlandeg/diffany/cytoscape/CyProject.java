package be.svlandeg.diffany.cytoscape;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.CyNetworkViewManager;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.internal.Services;
import be.svlandeg.diffany.semantics.*;



/**
 * The CyProject keeps the GUI settings for the current project and is able to generate a 
 * {@link Project} object to execute the algorithm.
 * 
 * @author Thomas Van Parys
 *
 */
public class CyProject{

	private Project project;
	
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
	private Set<CyNetworkPair> resultNetworks = new HashSet<CyNetworkPair>();;

	private EdgeOntology edgeOntology = new DefaultEdgeOntology();
	private NodeMapper nodeMapper = new DefaultNodeMapper();

	/**
	 * Default project name
	 */
	public static final String DEFAULT_PROJECT_NAME = "New Project";
	private String name = DEFAULT_PROJECT_NAME;
	
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
	 * Generates a {@link Project} from current information.
	 * 
	 * @return a newly constructed {@link Project}
	 * @throws InvalidProjectException is thrown when not all necessary parameters are there to construct the {@link Project}
	 */
	public Project getProject() throws InvalidProjectException{
		if (!canExecute()){
			throw new InvalidProjectException();
		}
		
		CyNetworkBridge bridge = new CyNetworkBridge();
		ReferenceNetwork refNet = bridge.getReferenceNetwork(this.getReferenceNetwork(), this.edgeOntology);
		
		Set<ConditionNetwork> condNets = new HashSet<ConditionNetwork>();
		for (CyNetwork cyCondNet : this.getConditionalNetworks()){
			ConditionNetwork condNet = bridge.getConditionNetwork(cyCondNet, this.edgeOntology);
			condNets.add(condNet);			
		}
		
		this.project = new Project(this.getName(), refNet, condNets, this.getEdgeOntology(), this.getNodeMapper());
		
		return project;
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
	 * Get the {@link EdgeOntology}.
	 * 
	 * @return the edge ontology defined for this project
	 */
	public EdgeOntology getEdgeOntology() {
		return edgeOntology;
	}

	/**
	 * Set the {@link EdgeOntology}.
	 * @param edgeOntology the edge ontology defined for this project
	 */
	public void setEdgeOntology(EdgeOntology edgeOntology) {
		this.edgeOntology = edgeOntology;
	}

	/**
	 * Get the {@link NodeMapper}.
	 * @return node mapper defined for this project
	 */
	public NodeMapper getNodeMapper() {
		return nodeMapper;
	}

	/**
	 * Set the {@link NodeMapper}.
	 * @param nodeMapper node mapper defined for this project
	 */
	public void setNodeMapper(NodeMapper nodeMapper) {
		this.nodeMapper = nodeMapper;
	}

	/**
	 * Get the project name
	 * @return the project name
	 */
	public String getName() {
		return name;
	}

	/**
	 * Set the project name
	 * @param name the project name
	 */
	public void setName(String name) {
		this.name = name;
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
