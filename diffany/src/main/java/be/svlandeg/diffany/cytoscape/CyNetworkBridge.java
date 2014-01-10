package be.svlandeg.diffany.cytoscape;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyIdentifiable;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyRow;
import org.cytoscape.model.CyTable;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CySubNetwork;

import be.svlandeg.diffany.concepts.Condition;
import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Node;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.semantics.EdgeOntology;

/**
 * Conversion class that takes the appropriate Cytoscape Factories and uses them
 * to create given concepts. The other way around, it enables existing {@link CyNetwork}s to
 * be converted to the project's {@link Network} concept. 
 * 
 * @author Thomas Van Parys
 *
 */
public class CyNetworkBridge {
	

	/**
	 * Column name of an extra field to be used to map nodes from different networks onto each other.
	 */
	public static final String NORMALIZED_NAME = "normalized_name.SUID";
	
	/**
	 * Column name for edge weight (double)
	 */
	public static String WEIGHT = "weight";
	/**
	 * Column name for negated edges (boolean)
	 */
	public static String NEGATED = "negated";
	
	/**
	 * Construct new bridge object
	 */
	public CyNetworkBridge(){
	}
	
	/**
	 * Creates a {@link CyNetwork} from a {@link Network} and adds it to Cytoscape.
	 * 
	 * @param network the {@link Network} to be added in Cytoscape
	 * @param collection the {@link CyRootNetwork} to link the new {@link Network} to
	 * @return the newly created and added {@link CyNetwork}
	 */
	public CyNetwork createCyNetwork(Network network, CyRootNetwork collection){
		CySubNetwork cyNetwork = collection.addSubNetwork();
		return createCyNetwork(network, cyNetwork);
	}
	
	/**
	 * Creates a {@link CyNetwork} from a {@link Network} and adds it to Cytoscape.
	 * 
	 * @param network the {@link Network} to be added in Cytoscape
	 * @param factory the {@link CyNetworkFactory} to create a fresh {@link CyNetwork}
	 * @return the newly created and added {@link CyNetwork}
	 */
	public CyNetwork createCyNetwork(Network network, CyNetworkFactory factory){
		return createCyNetwork(network, factory.createNetwork());
	}
	/**
	 * Converts a {@link Network} to a given {@link CyNetwork}
	 * 
	 * @param network input {@link Network} 
	 * @param cyNetwork the new network, created by the {@link CyNetworkFactory}
	 * @return {@link CyNetwork} that can be added to Cytoscape using a {@link CyNetworkManager}
	 */
	private CyNetwork createCyNetwork(Network network, CyNetwork cyNetwork){
		
		cyNetwork.getRow(cyNetwork).set(CyNetwork.NAME,network.getName());
		
		CyTable nodeTable = cyNetwork.getDefaultNodeTable();
		CyTable edgeTable = cyNetwork.getDefaultEdgeTable();

		if (nodeTable.getColumn(NORMALIZED_NAME)==null){
			nodeTable.createColumn(NORMALIZED_NAME, String.class, false);			
		}

		if (edgeTable.getColumn(WEIGHT)==null){
			edgeTable.createColumn(WEIGHT, Double.class, false);			
		}
		
		if (edgeTable.getColumn(NEGATED)==null){
			edgeTable.createColumn(NEGATED, Boolean.class, false);			
		}
		
		for (Node node: network.getNodes()){
			CyNode cyNode = cyNetwork.addNode();
			cyNetwork.getRow(cyNode).set(CyNetwork.NAME, node.getName());
			cyNetwork.getRow(cyNode).set(NORMALIZED_NAME, node.getName(true));
		}
		
		for (Edge edge : network.getEdges()){
			Node from = edge.getSource();
			Node to = edge.getTarget();
			
			CyTable table = cyNetwork.getDefaultNodeTable();
			final String pk = table.getPrimaryKey().getName();
			CyNode fromNode = nodeToCyNode(from, cyNetwork, table, pk);
			CyNode toNode = nodeToCyNode(to, cyNetwork, table, pk);
			
			CyEdge newEdge = cyNetwork.addEdge(fromNode, toNode, !edge.isSymmetrical());
			cyNetwork.getRow(newEdge).set(WEIGHT, edge.getWeight());
			cyNetwork.getRow(newEdge).set(CyEdge.INTERACTION, edge.getType());
			cyNetwork.getRow(newEdge).set(NEGATED, edge.isNegated());
			
			
		}
		return cyNetwork;
	}
	
	/**
	 * Takes a {@link Node} and converts it to a {@link CyNode}, using the normalized_name column.
	 * 
	 * @param node
	 * @param cyNetwork
	 * @param table
	 * @param pk
	 * @return
	 */
	private CyNode nodeToCyNode(Node node, CyNetwork cyNetwork, CyTable table, String pk){
		CyNode cyNode = null;
		Collection<CyRow> matches = table.getMatchingRows(NORMALIZED_NAME, node.getName(true));
		for (final CyRow row : matches){
			Long suid = row.get(pk, Long.class);
			cyNode = cyNetwork.getNode(suid);
		}
		return cyNode;
	}
	
	/**
	 * To distinguish between the two types of networks
	 * 
	 * @author Thomas Van Parys
	 *
	 */
	public enum NetworkType{
		REFERENCE, CONDITION;
	}
	
	/**
	 * Gets the {@link CyNetwork} from Cytoscape
	 * and creates a {@link Network} out of it.
	 * 
	 * @param network the "shared name" of the network in the network table.
	 * 
	 * 
	 * @return the equivalent {@link Network} object
	 */
	public ReferenceNetwork getReferenceNetwork(CyNetwork network, EdgeOntology edgeOntology){
		return (ReferenceNetwork)getNetwork(network, NetworkType.REFERENCE, edgeOntology);
	}
	
	/**
	 * Gets the {@link CyNetwork} from Cytoscape
	 * and creates a {@link Network} out of it.
	 * 
	 * @param network the "shared name" of the network in the network table.
	 * 
	 * 
	 * @return the equivalent {@link Network} object
	 */
	public ConditionNetwork getConditionNetwork(CyNetwork network, EdgeOntology edgeOntology){
		return (ConditionNetwork)getNetwork(network, NetworkType.CONDITION, edgeOntology);
	}
	
	/**
	 * Gets the {@link CyNetwork} from Cytoscape, based on the name,
	 * and creates a {@link Network} out of it.
	 * 
	 * @param cyNetwork the "name" of the network in the network table.
	 * @param type a switch to determine with role this network plays in the {@link Project}
	 * @param edgeOntology will determine which edge types are considered symmetrical
	 * 
	 * @return the equivalent {@link Network} object
	 */
	private Network getNetwork(CyNetwork cyNetwork, NetworkType type, EdgeOntology edgeOntology){
		
		Network network = null;
		String netName = this.getName(cyNetwork, cyNetwork);
		switch(type){
		case REFERENCE: 
			network = new ReferenceNetwork(netName);
			break;
		case CONDITION:
			//TODO get conditions from gui
			Set<Condition> conditions = new HashSet<Condition>();
			conditions.add(new Condition("temp_condition"));
			network = new ConditionNetwork(netName, conditions);
			break;
		}
				
		List<CyNode> cyNodes = cyNetwork.getNodeList();
		Map<Long, Node> nodeMap = new HashMap<Long, Node>();
		for (CyNode cyNode : cyNodes){
			String nodeName = this.getName(cyNetwork, cyNode);
			Node node = new Node(nodeName);
			nodeMap.put(cyNode.getSUID(), node);
		}
		
		Set<Edge> edgeSet = new HashSet<Edge>();
		List<CyEdge> cyEdges = cyNetwork.getEdgeList();
		for (CyEdge cyEdge : cyEdges){
			Long nodeFromID = cyEdge.getSource().getSUID();
			Long nodeToID = cyEdge.getTarget().getSUID();
			Node fromNode = nodeMap.get(nodeFromID);
			Node toNode = nodeMap.get(nodeToID);
			
			String interaction = this.getInteraction(cyNetwork, cyEdge);
			
			boolean isSymmetrical = edgeOntology.isSymmetricalSourceType(interaction);
			
			Edge edge = new Edge(interaction, fromNode, toNode, isSymmetrical);
			if (this.getWeight(cyNetwork, cyEdge) !=null){
				edge.setWeight(this.getWeight(cyNetwork, cyEdge));				
			}
			edge.makeNegated(this.isNegated(cyNetwork, cyEdge));
			edgeSet.add(edge);
		}
		
		network.setNodesAndEdges(new HashSet<Node>(nodeMap.values()), edgeSet);
		
		return network;
	}
	
	/**
	 * Get the value of the NEGATED column
	 * 
	 * @param cyNetwork
	 * @param cyEdge
	 * @return true if the NEGATED column contains "true"
	 */
	private boolean isNegated(CyNetwork cyNetwork, CyEdge cyEdge) {
		Boolean negated = cyNetwork.getRow(cyEdge).get(NEGATED, Boolean.class);
		return negated !=null && negated;
	}

	/**
	 * Get the value of the NAME column.
	 * 
	 * @param cyNetwork the network containing the node
	 * @param cyObject the object 
	 * @return a string representing the object's name
	 */
	private String getName(CyNetwork cyNetwork, CyIdentifiable cyObject) {
		return cyNetwork.getRow(cyObject).get(CyNetwork.NAME, String.class);
	}
	/**
	 * Get the value of the INTERACTION column
	 * 
	 * @param cyNetwork
	 * @param cyEdge
	 * @return
	 */
	private String getInteraction(CyNetwork cyNetwork, CyEdge cyEdge){
		return cyNetwork.getRow(cyEdge).get(CyEdge.INTERACTION, String.class);
	}
	/**
	 * Get the value of the WEIGHT column
	 * 
	 * @param cyNetwork
	 * @param cyEdge
	 * @return the value of the weight column, or null if the column doesn't exist.
	 */
	private Double getWeight(CyNetwork cyNetwork, CyEdge cyEdge){
		return cyNetwork.getRow(cyEdge).get(this.WEIGHT, Double.class);
	}
}
