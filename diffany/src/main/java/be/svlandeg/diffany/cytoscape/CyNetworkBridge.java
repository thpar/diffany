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

import be.svlandeg.diffany.concepts.Condition;
import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.EdgeDefinition;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Node;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.internal.Services;

/**
 * Conversion class that takes the appropriate Cytoscape Factories and uses them
 * to create given concepts. The other way around, it enables existing {@link CyNetwork}s to
 * be converted to the project's {@link Network} concept. 
 * 
 * @author thpar
 *
 */
public class CyNetworkBridge {
	

	/**
	 * Column name of an extra field to be used to map nodes from different networks onto each other.
	 */
	private final String NORMALIZED_NAME = "normalized_name";

	private Services services;

	private Model model;

	public CyNetworkBridge(Model model){
		this.services = model.getServices();
		this.model = model;
	}

	/**
	 * Converts a {@link Network} to a {@link CyNetwork}
	 * 
	 * @param network input {@link Network} 
	 * @return {@link CyNetwork} that can be added to Cytoscape using a {@link CyNetworkManager}
	 */
	public CyNetwork createCyNetwork(Network network){
		CyNetworkFactory cyNetworkFactory = services.getCyNetworkFactory();
		CyNetwork cyNetwork = cyNetworkFactory.createNetwork();
		cyNetwork.getRow(cyNetwork).set(CyNetwork.NAME,network.getName());
		cyNetwork.getDefaultNodeTable().createColumn(NORMALIZED_NAME, String.class, false);
		
		
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
			
			cyNetwork.addEdge(fromNode, toNode, !edge.isSymmetrical());
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
	
	enum NetworkType{
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
	public ReferenceNetwork getReferenceNetwork(CyNetwork network){
		return (ReferenceNetwork)getNetwork(network, NetworkType.REFERENCE);
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
	public ConditionNetwork getConditionNetwork(CyNetwork network){
		return (ConditionNetwork)getNetwork(network, NetworkType.CONDITION);
	}
	
	/**
	 * Gets the {@link CyNetwork} from Cytoscape, based on the name,
	 * and creates a {@link Network} out of it.
	 * 
	 * @param cyNetwork the "name" of the network in the network table.
	 * 
	 * 
	 * @return the equivalent {@link Network} object
	 */
	private Network getNetwork(CyNetwork cyNetwork, NetworkType type){
		
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
			
			EdgeDefinition def = new EdgeDefinition();
			def.setType(this.getInteraction(cyNetwork, cyEdge));
			Edge edge = new Edge(fromNode, toNode, def);
			edgeSet.add(edge);
		}
		
		network.setNodesAndEdges(new HashSet<Node>(nodeMap.values()), edgeSet);
		
		return network;
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
	
}
