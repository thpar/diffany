package be.svlandeg.diffany.cytoscape;

import java.util.HashSet;
import java.util.Set;

import org.cytoscape.model.CyNetwork;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;



/**
 * The CyProject keeps the GUI settings for the current project and is able to generate a 
 * {@link Project} object to execute the algorithm.
 * 
 * @author thpar
 *
 */
public class CyProject{

	private Project project;
	
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

	public static final String DEFAULT_PROJECT_NAME = "New Project";
	private String name = DEFAULT_PROJECT_NAME;
	
	public CyProject(){
	}
	
	
	public CyNetwork getReferenceNetwork() {
		return referenceNetwork;
	}



	public void setReferenceNetwork(CyNetwork referenceNetwork) {
		this.referenceNetwork = referenceNetwork;
	}

	public void removeConditionalNetwork(CyNetwork cyNetwork){
		this.conditionalNetworks.remove(cyNetwork);
	}


	public Set<CyNetwork> getConditionalNetworks() {
		return conditionalNetworks;
	}



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
		ReferenceNetwork refNet = bridge.getReferenceNetwork(this.getReferenceNetwork());
		
		Set<ConditionNetwork> condNets = new HashSet<ConditionNetwork>();
		for (CyNetwork cyCondNet : this.getConditionalNetworks()){
			ConditionNetwork condNet = bridge.getConditionNetwork(cyCondNet);
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
		//TODO actual check 
		return true;
	}


	public EdgeOntology getEdgeOntology() {
		return edgeOntology;
	}


	public void setEdgeOntology(EdgeOntology edgeOntology) {
		this.edgeOntology = edgeOntology;
	}


	public NodeMapper getNodeMapper() {
		return nodeMapper;
	}


	public void setNodeMapper(NodeMapper nodeMapper) {
		this.nodeMapper = nodeMapper;
	}


	public String getName() {
		return name;
	}


	public void setName(String name) {
		this.name = name;
	}


	public void addResultPair(CyNetwork diffNet, CyNetwork overlapNet) {
		this.resultNetworks.add(new CyNetworkPair(diffNet, overlapNet));
	}
	
	public Set<CyNetworkPair> getResultNetworks() {
		return resultNetworks;
	}
		

}
