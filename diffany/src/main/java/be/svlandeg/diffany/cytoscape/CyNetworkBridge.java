package be.svlandeg.diffany.cytoscape;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.cytoscape.model.CyColumn;
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
import org.cytoscape.view.layout.CyLayoutAlgorithm;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskManager;

import be.svlandeg.diffany.core.networks.Condition;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.EdgeDefinition;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.semantics.EdgeOntology;
import be.svlandeg.diffany.cytoscape.internal.Services;

/**
 * Static conversion class that takes the appropriate Cytoscape Factories and uses them
 * to create given concepts. The other way around, it enables existing {@link CyNetwork}s to
 * be converted to the project's {@link Network} concept. 
 * 
 * @author Thomas Van Parys
 *
 */
public class CyNetworkBridge {
	
	
	/**
	 * Column name of an extra field to be used as preferred node label (probably the symbolic gene name)
	 */
	public static final String SYMBOLIC_NAME = "official_symbol";
	
	/**
	 * Column name for edge weight (double)
	 */
	public static final String WEIGHT = "weight";
	/**
	 * Column name for negated edges (boolean)
	 */
	public static final String NEGATED = "negated";
	
	/**
	 * Unique ID for the next network
	 */
	private static int nextID = 0;
	
	/**
	 * {@link CyNode} attributes that should not be transferred to the {@link Node} attributes.
	 */
	static final Set<String> defaultNodeAttrs = new HashSet<String>();
	static {
		defaultNodeAttrs.add("SUID");
		defaultNodeAttrs.add("selected");
		defaultNodeAttrs.add("shared name");
		defaultNodeAttrs.add(SYMBOLIC_NAME);
		defaultNodeAttrs.add("name");
	}
	
	/**
	 * Creates a {@link CyNetwork} from a {@link Network} and adds it to Cytoscape.
	 * 
	 * @param network the {@link Network} to be added in Cytoscape
	 * @param collection the {@link CyRootNetwork} to put new {@link CyNetwork} in
	 * @return the newly created and added {@link CyNetwork}
	 */
	public static CyNetwork createCyNetwork(Network network, CyRootNetwork collection){
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
	public static CyNetwork createCyNetwork(Network network, CyNetworkFactory factory){
		return createCyNetwork(network, factory.createNetwork());
	}
	
	/**
	 * Converts a {@link Network} to a given {@link CyNetwork}
	 * 
	 * @param network input {@link Network} 
	 * @param cyNetwork the new network, created by a {@link CyNetworkFactory} or {@link CyRootNetwork}
	 * @return {@link CyNetwork} that can be added to Cytoscape using a {@link CyNetworkManager}
	 */
	private static CyNetwork createCyNetwork(Network network, CyNetwork cyNetwork){
		
		cyNetwork.getRow(cyNetwork).set(CyNetwork.NAME,network.getName());
		
		CyTable nodeTable = cyNetwork.getTable(CyNode.class, CyNetwork.LOCAL_ATTRS);
		CyTable edgeTable = cyNetwork.getTable(CyEdge.class, CyNetwork.LOCAL_ATTRS);
		
		if (nodeTable.getColumn(SYMBOLIC_NAME)==null){
			nodeTable.createColumn(SYMBOLIC_NAME, String.class, false);			
		}

		if (edgeTable.getColumn(WEIGHT)==null){
			edgeTable.createColumn(WEIGHT, Double.class, false);			
		}
		
		if (edgeTable.getColumn(NEGATED)==null){
			edgeTable.createColumn(NEGATED, Boolean.class, false);			
		}
		
		//add Network specific columns for all extra node attributes
		for (String attr : network.getAllNodeAttributes()){
			if (nodeTable.getColumn(attr)==null){
				nodeTable.createColumn(attr, String.class, false);			
			}
		}
		
		for (Node node: network.getNodes()){
			CyNode cyNode = cyNetwork.addNode();
			cyNetwork.getRow(cyNode).set(CyNetwork.NAME, node.getID());
			cyNetwork.getRow(cyNode).set(SYMBOLIC_NAME, node.getDisplayName());
			
			for (String attr : network.getAllNodeAttributes()){
				cyNetwork.getRow(cyNode).set(attr, node.getAttribute(attr));
			}
		}
		
		for (Edge edge : network.getEdges()){
			Node from = edge.getSource();
			Node to = edge.getTarget();
			
			CyNode fromNode = findCyNodeByNode(from, cyNetwork);
			CyNode toNode = findCyNodeByNode(to, cyNetwork);
			
			CyEdge newEdge = cyNetwork.addEdge(fromNode, toNode, !edge.isSymmetrical());
			cyNetwork.getRow(newEdge).set(WEIGHT, edge.getWeight());
			cyNetwork.getRow(newEdge).set(CyEdge.INTERACTION, edge.getType());
			cyNetwork.getRow(newEdge).set(NEGATED, edge.isNegated());
			
			
		}
		return cyNetwork;
	}
	
	/**
	 * Takes a {@link Node} and looks it up in the given {@link CyNetwork}, using the unique nodeID.
	 * If no nodeID is returned, the normalized_name column will be used against a normalized version of the node name.
	 * 
	 * @param node {@link Node} to be converted
	 * @param cyNetwork the {@link CyNetwork} to search in
	 * @return the {@link CyNode} with the same normalized name as the {@link Node}.
	 * Returns null if no match was found.
	 */
	private static CyNode findCyNodeByNode(Node node, CyNetwork cyNetwork){
		CyTable table = cyNetwork.getDefaultNodeTable();
		final String pk = table.getPrimaryKey().getName();
		
		CyNode cyNode = null;
		
		Collection<CyRow> matches = table.getMatchingRows(CyNetwork.NAME, node.getID());
		
		for (final CyRow row : matches){
			Long suid = row.get(pk, Long.class);
			cyNode = cyNetwork.getNode(suid);
		}
		return cyNode;
	}
	
	/**
	 * To distinguish between two types of source networks
	 * 
	 * @author Thomas Van Parys
	 *
	 */
	public enum NetworkType{
		REFERENCE, CONDITION, GENERIC;
	}
	
	/**
	 * Gets the {@link CyNetwork} from Cytoscape and creates a {@link ReferenceNetwork} out of it.
	 * 
	 * @param network the "shared name" of the network in the network table.
	 * @param edgeOntology will determine which edge types are considered symmetrical
	 * 
	 * @return the equivalent {@link Network} object
	 */
	public static ReferenceNetwork getReferenceNetwork(CyNetwork network, EdgeOntology edgeOntology){
		return (ReferenceNetwork)getNetwork(network, NetworkType.REFERENCE, edgeOntology);
	}
	
	/**
	 * Gets the {@link CyNetwork} from Cytoscape and creates a {@link ConditionNetwork} out of it.
	 * 
	 * @param network the "shared name" of the network in the network table.
	 * @param edgeOntology will determine which edge types are considered symmetrical
	 * 
	 * @return the equivalent {@link Network} object
	 */
	public static ConditionNetwork getConditionNetwork(CyNetwork network, EdgeOntology edgeOntology){
		return (ConditionNetwork)getNetwork(network, NetworkType.CONDITION, edgeOntology);
	}
	
	/**
	 * Gets the {@link CyNetwork} from Cytoscape and creates a {@link InputNetwork} out of it.
	 * 
	 * @param network the "shared name" of the network in the network table.
	 * @param edgeOntology will determine which edge types are considered symmetrical
	 * 
	 * @return the equivalent {@link Network} object
	 */
	public static InputNetwork getInputNetwork(CyNetwork network, EdgeOntology edgeOntology){
		return (InputNetwork)getNetwork(network, NetworkType.GENERIC, edgeOntology);
	}
	
	/**
	 * Gets the {@link CyNetwork} from Cytoscape, based on the name,
	 * and creates a {@link Network} out of it.
	 * 
	 * @param cyNetwork the "name" of the network in the network table.
	 * @param type a switch to determine which role this network plays in the {@link Project}
	 * @param edgeOntology will determine which edge types are considered symmetrical
	 * 
	 * @return the equivalent {@link Network} object
	 */
	private static Network getNetwork(CyNetwork cyNetwork, NetworkType type, EdgeOntology edgeOntology){
	
		Network network = null;
		String netName = getName(cyNetwork, cyNetwork);
		
		//extract extra node columns
		Set<String> nodeAttrs = new HashSet<String>();
		CyTable table = cyNetwork.getDefaultNodeTable();
		Collection<CyColumn> cyCols = table.getColumns();
		for (CyColumn cyCol : cyCols){
			String colName = cyCol.getName();
			if (!defaultNodeAttrs.contains(colName)){
				nodeAttrs.add(colName);				
			}
		}
		
		switch(type){
		case REFERENCE: 
			network = new ReferenceNetwork(netName, getNextID(), nodeAttrs);
			break;
		case CONDITION:
			//TODO v3.0: get conditions from gui
			Set<Condition> conditions = new HashSet<Condition>();
			conditions.add(new Condition("undefined_condition"));
			network = new ConditionNetwork(netName, getNextID(), nodeAttrs, conditions);
			break;
		case GENERIC:
			network = new InputNetwork(netName, getNextID(), nodeAttrs);
			break;
		}
				
		List<CyNode> cyNodes = cyNetwork.getNodeList();
		Map<Long, Node> nodeMap = new HashMap<Long, Node>();
		for (CyNode cyNode : cyNodes){
			String nodeID = getName(cyNetwork, cyNode);
			String symbolicName = getSymbolicName(cyNetwork, cyNode);
			if (symbolicName==null || symbolicName.isEmpty()){
				symbolicName = nodeID;
			}
			Node node = new Node(nodeID, symbolicName);
			
			//transfer node attributes
			for (String attr : nodeAttrs){
				node.setAttribute(attr, getAttribute(cyNetwork, cyNode, attr));
			}
			
			nodeMap.put(cyNode.getSUID(), node);	
		}
		
		Set<Edge> edgeSet = new HashSet<Edge>();
		List<CyEdge> cyEdges = cyNetwork.getEdgeList();
		for (CyEdge cyEdge : cyEdges){
			Long nodeFromID = cyEdge.getSource().getSUID();
			Long nodeToID = cyEdge.getTarget().getSUID();
			Node fromNode = nodeMap.get(nodeFromID);
			Node toNode = nodeMap.get(nodeToID);
			
			String interaction = getInteraction(cyNetwork, cyEdge);
			
			boolean isSymmetrical = edgeOntology.isSymmetricalSourceType(interaction);
			boolean isNegated = isNegated(cyNetwork, cyEdge);
			Double weight = getWeight(cyNetwork, cyEdge); 
			if (weight == null){
				weight = EdgeDefinition.DEFAULT_WEIGHT;				
			}
			
			Edge edge = new Edge(interaction, fromNode, toNode, isSymmetrical, weight, isNegated);
			
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
	private static boolean isNegated(CyNetwork cyNetwork, CyEdge cyEdge) {
		Boolean negated = cyNetwork.getRow(cyEdge).get(NEGATED, Boolean.class);
		return negated !=null && negated;
	}

	/**
	 * Get the value of the NAME column. Works for all CyObjects with a name.
	 * 
	 * @param cyNetwork the network containing the node
	 * @param cyObject the object 
	 * @return a string representing the object's name
	 */
	private static String getName(CyNetwork cyNetwork, CyIdentifiable cyObject) {
		return cyNetwork.getRow(cyObject).get(CyNetwork.NAME, String.class);
	}
	
	/**
	 * Get the symbolic name of the node
	 * 
	 * @param cyNetwork the network containing the node
	 * @param cyNode the node
	 * @return a string representing the symbolic name of the node (used for labeling)
	 * 
	 */
	private static String getSymbolicName(CyNetwork cyNetwork, CyNode cyNode) {
		return cyNetwork.getRow(cyNode).get(SYMBOLIC_NAME, String.class);
	}
	
	/**
	 * Get the value of an attribute column.
	 * Values will always be Strings, even if the column originally contains int or bool
	 * 
	 * @param cyNetwork the network containing the node
	 * @param cyNode the node
	 * @param columnName name of the attribute column
	 * @return value of the column for this specific {@link CyNode}
	 */
	private static String getAttribute(CyNetwork cyNetwork, CyNode cyNode, String columnName) {
		String attr = cyNetwork.getRow(cyNode).get(columnName, String.class);
		if (attr == null){
			return new String();
		} else {
			return attr.toString();
		}
		
	}
	
	/**
	 * Get the value of the INTERACTION column for a given edge.
	 * 
	 * @param cyNetwork
	 * @param cyEdge
	 * @return
	 */
	private static String getInteraction(CyNetwork cyNetwork, CyEdge cyEdge){
		return cyNetwork.getRow(cyEdge).get(CyEdge.INTERACTION, String.class);
	}
	/**
	 * Get the value of the WEIGHT column for a given edge.
	 * 
	 * @param cyNetwork
	 * @param cyEdge
	 * @return the value of the weight column, or null if the column doesn't exist.
	 */
	private static Double getWeight(CyNetwork cyNetwork, CyEdge cyEdge){
		return cyNetwork.getRow(cyEdge).get(WEIGHT, Double.class);
	}
	
	
	/**
	 * Converts a {@link Network} to a {@link CyNetwork} and registers it
	 * with Cytoscape as a {@link CySubNetwork} of the given Collection.
	 * 
	 * @param network the {@link Network} to be added
	 * @param collection the {@link CyRootNetwork} the new network should belong to.
	 * @param services the {@link Services} from the context
	 * @return the created and added {@link CyNetwork}
	 */
	public static CyNetwork addCyNetwork(Network network, CyRootNetwork collection, Services services){
		CyNetwork cyNet = createCyNetwork(network, collection);
		nextID++; //keep track of ID's handed out when adding new networks
		
		services.getCyNetworkManager().addNetwork(cyNet);
		
		CyNetworkView cyView = services.getCyNetworkViewFactory().createNetworkView(cyNet);
		services.getCyNetworkViewManager().addNetworkView(cyView);
		
		CyLayoutAlgorithm layout = services.getCyLayoutAlgorithmManager().getLayout("force-directed");
		TaskIterator it = layout.createTaskIterator(cyView, layout.createLayoutContext(), CyLayoutAlgorithm.ALL_NODE_VIEWS, null);
		TaskManager<?, ?> tm = services.getTaskManager();
		tm.execute(it);
		
		return cyNet;
	}
	
	/**
	 * 
	 * @return an ID unique for this session
	 */
	public static int getNextID(){
		return nextID++;
	}
	
}
