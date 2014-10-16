package be.svlandeg.diffany.core.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.activation.UnsupportedDataTypeException;

import be.svlandeg.diffany.core.networks.Condition;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.semantics.NodeMapper;

/**
 * This class allows reading or writing a {@link Network} from File.
 * 
 * @author Sofie Van Landeghem
 */
public class NetworkIO
{

	private static String default_conditions_file = "conditions.txt";
	private static String default_edge_file = "edges.txt";
	private static String default_node_file = "nodes.txt";
	private static String default_definition_file = "network.txt";

	private static String name_field = "Name";
	private static String ID_field = "ID";
	private static String type_field = "Type";
	private static String attributes_field = "Attributes";

	/**
	 * Write the string representation of the network: All edges to one File, and all nodes to another, as specified by the parameters.
	 * Additionally, the network name, type and node attributes are written to a defitionsfile.
	 * 
	 * @param network the {@link Network} that needs to be written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param edgesFile the output file in which the edges will be written
	 * @param nodesFile the output file in which the nodes will be written
	 * @param definitionFile the output file in which the network definition will be written
	 * @param writeHeaders whether or not to write the headers of the nodes and edges files
	 * @param allowVirtualEdges if true, write virtual edges in the edges file. If false, condense the information into node attributes
	 * @see EdgeIO#writeToTab
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	protected static void writeNetworkToFiles(Network network, NodeMapper nm, File edgesFile, File nodesFile, File definitionFile, boolean writeHeaders, boolean allowVirtualEdges) throws IOException
	{
		// EDGES
		edgesFile.getParentFile().mkdirs();
		BufferedWriter edgeWriter = new BufferedWriter(new FileWriter(edgesFile));

		if (writeHeaders)
		{
			edgeWriter.append(EdgeIO.getHeader());
			edgeWriter.newLine();
			edgeWriter.flush();
		}

		Set<Edge> normalEdges = network.getEdgesByVirtualState(false);
		Set<Edge> virtualEdges = network.getEdgesByVirtualState(true);

		for (Edge e : normalEdges)
		{
			edgeWriter.append(EdgeIO.writeToTab(e));
			edgeWriter.newLine();
			edgeWriter.flush();
		}

		if (allowVirtualEdges)
		{
			for (Edge e : virtualEdges)
			{
				edgeWriter.append(EdgeIO.writeToTab(e));
				edgeWriter.newLine();
				edgeWriter.flush();
			}
		}
		edgeWriter.flush();
		edgeWriter.close();

		// NODES
		nodesFile.getParentFile().mkdirs();
		BufferedWriter nodeWriter = new BufferedWriter(new FileWriter(nodesFile));

		Set<Node> unduplicatedNodes = new HashSet<Node>();
		for (Node n : network.getNodes())
		{
			if (!nm.isContained(n, unduplicatedNodes))
			{
				// exclude virtual nodes when virtual edges should not be written
				if (allowVirtualEdges || !n.isVirtual())
				{
					unduplicatedNodes.add(n);
				}
			}
		}

		SortedSet<String> nodeAttributes = new TreeSet<String>();
		nodeAttributes.addAll(network.getAllNodeAttributes());
		
		if (writeHeaders)
		{
			nodeWriter.append(NodeIO.getHeader(nodeAttributes));
			nodeWriter.newLine();
			nodeWriter.flush();
		}

		SortedMap<String, Node> sortedNodes = new TreeMap<String, Node>();
		for (Node n : unduplicatedNodes)
		{
			sortedNodes.put(n.getID(), n);
		}
		for (String ID : sortedNodes.keySet())
		{
			Node n = sortedNodes.get(ID);
			String tabRep = NodeIO.writeToTab(n, nodeAttributes);
			nodeWriter.append(tabRep);
			nodeWriter.newLine();
			nodeWriter.flush();
		}
		nodeWriter.flush();
		nodeWriter.close();

		// DEFINITION: NAME, ID and TYPE (CLASS)
		// TODO create DefinitionIO to read and write the network definition
		definitionFile.getParentFile().mkdirs();
		BufferedWriter defWriter = new BufferedWriter(new FileWriter(definitionFile));

		defWriter.append(ID_field + "\t" + network.getID());
		defWriter.newLine();

		defWriter.append(name_field + "\t" + network.getName());
		defWriter.newLine();

		String networkClass = network.getClass().getSimpleName();
		defWriter.append(type_field + "\t" + networkClass);
		defWriter.newLine();

		String attributes = "";
		for (String attribute : network.getAllNodeAttributes())
		{
			attributes += attribute + ";";
		}
		defWriter.append(attributes_field + "\t" + attributes);
		defWriter.newLine();

		defWriter.flush();
		defWriter.close();
	}

	/**
	 * Write the string representation of the conditions of a condition-specific network: call writeNetworkToFiles to write the edges edges information, nodes and network definition first.
	 * Subsequently, call this method to write the conditions (ontology terms).
	 * 
	 * @param network the {@link ConditionNetwork} that needs to be written
	 * @param conditionsFile the file to which the input network should be written
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	public static void writeConditionsToFiles(ConditionNetwork network, File conditionsFile) throws IOException
	{
		Set<Condition> conditions = network.getConditions();

		conditionsFile.getParentFile().mkdirs();
		BufferedWriter writer = new BufferedWriter(new FileWriter(conditionsFile));
		for (Condition c : conditions)
		{
			writer.append(c.getDescription());
			for (String ont : c.getOntologyTerms())
			{
				writer.append("\t" + ont);
			}
			writer.newLine();
			writer.flush();
		}

		writer.flush();
		writer.close();
	}

	/**
	 * Write the string representation of the network to (subfiles of) a directory: edges, nodes and definition all go in a separate file.
	 * If the network is a ConditionNetwork, the conditions will be written to a fourth file.
	 * 
	 * @param network the {@link Network} that needs to be written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param dir the output dir in which the tab files will be written
	 * @param writeHeaders whether or not to write the headers of the nodes and edges files
	 * @param allowVirtualEdges if true, write virtual edges in the edges file. If false, condense the information into node attributes
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	public static void writeNetworkToDir(Network network, NodeMapper nm, File dir, boolean writeHeaders, boolean allowVirtualEdges) throws IOException
	{
		File edgeFile = new File(dir.getAbsolutePath() + "/" + default_edge_file);
		File nodeFile = new File(dir.getAbsolutePath() + "/" + default_node_file);
		File definitionFile = new File(dir.getAbsolutePath() + "/" + default_definition_file);
		writeNetworkToFiles(network, nm, edgeFile, nodeFile, definitionFile, writeHeaders, allowVirtualEdges);

		if (network instanceof ConditionNetwork)
		{
			File default_conditions_File = new File(dir.getAbsolutePath() + "/" + default_conditions_file);
			writeConditionsToFiles((ConditionNetwork) network, default_conditions_File);
		}
	}

	/**
	 * Read a network from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a Network representation of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	private static Network readInputNetworkFromDir(File dir, NodeMapper nm, boolean skipHeader) throws IOException
	{
		File edgeFile = new File(dir.getAbsolutePath() + "/" + default_edge_file);
		File nodeFile = new File(dir.getAbsolutePath() + "/" + default_node_file);
		File definitionFile = new File(dir.getAbsolutePath() + "/" + default_definition_file);
		File conditionsFile = new File(dir.getAbsolutePath() + "/" + default_conditions_file);
		
		int ID = readIDFromFile(definitionFile);
		String name = readNameFromFile(definitionFile);
		String type = readTypeFromFile(definitionFile);
		List<String> listedAttributes = readAttributesFromFile(definitionFile);
		Set<String> attributes = new HashSet<String>(listedAttributes);

		Set<Node> nodes = readNodesFromFile(nodeFile, nm, skipHeader, listedAttributes);
		Set<Edge> edges = readEdgesFromFile(edgeFile, getMappedNodes(nodes), skipHeader);

		if (type.equals("ReferenceNetwork"))
		{
			ReferenceNetwork r = new ReferenceNetwork(name, ID, attributes, nm);
			r.setNodesAndEdges(nodes, edges);
			return r;
		}

		else if (type.equals("ConditionNetwork"))
		{
			Set<Condition> conditions = readConditionsFromFile(conditionsFile);
			ConditionNetwork c = new ConditionNetwork(name, ID, attributes, conditions, nm);
			c.setNodesAndEdges(nodes, edges);
			return c;
		}
		else if (type.equals("InputNetwork"))
		{
			InputNetwork c = new InputNetwork(name, ID, attributes, nm);
			c.setNodesAndEdges(nodes, edges);
			return c;
		}

		throw new UnsupportedDataTypeException("Encountered unknown input network type: " + type);
	}

	/**
	 * Read a {@link ReferenceNetwork} from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a ReferenceNetwork representation of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static ReferenceNetwork readReferenceNetworkFromDir(File dir, NodeMapper nm, boolean skipHeader) throws IOException
	{
		return (ReferenceNetwork) readInputNetworkFromDir(dir, nm, skipHeader);
	}

	/**
	 * Read a {@link ConditionNetwork} from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a ConditionNetwork representation of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static ConditionNetwork readConditionNetworkFromDir(File dir, NodeMapper nm, boolean skipHeader) throws IOException
	{
		return (ConditionNetwork) readInputNetworkFromDir(dir, nm, skipHeader);
	}

	/**
	 * Read an {@link InputNetwork} from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return an InputNetwork representation of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static InputNetwork readGenericInputNetworkFromDir(File dir, NodeMapper nm, boolean skipHeader) throws IOException
	{
		return (InputNetwork) readInputNetworkFromDir(dir, nm, skipHeader);
	}

	/**
	 * Read a set of {@link InputNetwork} from a directory, one network per subdirectory.
	 * 
	 * @param dir the output dir in which the subdirectories contain previously written tab files defining the different networks
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a set of InputNetwork representations of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static Set<InputNetwork> readGenericInputNetworksFromSubdirs(File dir, NodeMapper nm, boolean skipHeader) throws IOException
	{
		Set<InputNetwork> networks = new HashSet<InputNetwork>();
		for (File f : dir.listFiles())
		{
			if (f.isDirectory())
			{
				InputNetwork net = readGenericInputNetworkFromDir(f, nm, skipHeader);
				networks.add(net);
			}
		}
		return networks;
	}

	/**
	 * Read a network from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param reference the ReferenceNetwork linked to the read output network
	 * @param condNetworks the set of condition-specific networks linked to the read output network
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a Network representation of the nodes and edges in the files.
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	private static Network readOutputNetworkFromDir(File dir, NodeMapper nm, ReferenceNetwork reference, Set<ConditionNetwork> condNetworks, boolean skipHeader) throws IOException
	{
		File edgeFile = new File(dir.getAbsolutePath() + "/" + default_edge_file);
		File nodeFile = new File(dir.getAbsolutePath() + "/" + default_node_file);
		File definitionFile = new File(dir.getAbsolutePath() + "/" + default_definition_file);
		
		int ID = readIDFromFile(definitionFile);
		String name = readNameFromFile(definitionFile);
		String type = readTypeFromFile(definitionFile);
		List<String> listedAttributes = readAttributesFromFile(definitionFile);

		Set<Node> nodes = readNodesFromFile(nodeFile, nm, skipHeader, listedAttributes);
		Set<Edge> edges = readEdgesFromFile(edgeFile, getMappedNodes(nodes), skipHeader);


		if (type.equals("DifferentialNetwork"))
		{
			DifferentialNetwork d = new DifferentialNetwork(name, ID, reference, condNetworks, nm);
			d.setNodesAndEdges(nodes, edges);
			return d;
		}
		else if (type.equals("ConsensusNetwork"))
		{
			Set<Network> allNetworks = new HashSet<Network>();
			allNetworks.add(reference);
			allNetworks.addAll(condNetworks);
			ConsensusNetwork o = new ConsensusNetwork(name, ID, allNetworks, nm);
			o.setNodesAndEdges(nodes, edges);
			return o;
		}

		throw new UnsupportedDataTypeException("Encountered unknown output network type: " + type);
	}

	/**
	 * Retrieve a view on this set of nodes which is mapped by their unique IDs
	 * 
	 * @param nodes the original set of nodes
	 * @return the same set of nodes, mapped by their unique IDs
	 */
	private static Map<String, Node> getMappedNodes(Set<Node> nodes)
	{
		Map<String, Node> mappedNodes = new HashMap<String, Node>();
		for (Node n : nodes)
		{
			mappedNodes.put(n.getID(), n);
		}
		return mappedNodes;
	}

	/**
	 * Read a {@link DifferentialNetwork} from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param reference the ReferenceNetwork linked to this differential network
	 * @param condNetworks the set of condition-specific networks linked to this differential network
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a DifferentialNetwork representation of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static DifferentialNetwork readDifferentialNetworkFromDir(File dir, NodeMapper nm, ReferenceNetwork reference, Set<ConditionNetwork> condNetworks, boolean skipHeader) throws IOException
	{
		return (DifferentialNetwork) readOutputNetworkFromDir(dir, nm, reference, condNetworks, skipHeader);
	}

	/**
	 * Read a {@link ConsensusNetwork} from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param reference the ReferenceNetwork linked to this consensus network
	 * @param condNetworks the set of condition-specific networks linked to this consensus network
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a ConsensusNetwork representation of the nodes and edges in the files.
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static ConsensusNetwork readConsensusNetworkFromDir(File dir, NodeMapper nm, ReferenceNetwork reference, Set<ConditionNetwork> condNetworks, boolean skipHeader) throws IOException
	{
		return (ConsensusNetwork) readOutputNetworkFromDir(dir, nm, reference, condNetworks, skipHeader);
	}

	/**
	 * Read all conditions from a file containing one tab-delimited condition per line.
	 * 
	 * @param conditionsFile the file containing the condition data
	 * @return the set of conditions read from the file, or an empty set if no conditions were found
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static Set<Condition> readConditionsFromFile(File conditionsFile) throws IOException
	{
		Set<Condition> conditions = new HashSet<Condition>();

		BufferedReader reader = new BufferedReader(new FileReader(conditionsFile));
		String line = reader.readLine();
		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			String descr = stok.nextToken();
			Set<String> ontologies = new HashSet<String>();
			while (stok.hasMoreTokens())
			{
				ontologies.add(stok.nextToken());
			}
			Condition c = new Condition(descr, ontologies);
			conditions.add(c);

			line = reader.readLine();
		}

		reader.close();

		return conditions;
	}

	/**
	 * Read all edges from a file containing one tab-delimited edge per line. As input, a set of nodes should be given, mapped by their unique IDs.
	 * 
	 * @param edgesFile the file containing the edge data
	 * @param nodes the nodes relevant to the edges that will be read
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return the set of edges read from the file, or an empty set if no edges were found
	 * @see EdgeIO#readFromTab
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static Set<Edge> readEdgesFromFile(File edgesFile, Map<String, Node> nodes, boolean skipHeader) throws IOException
	{
		Set<Edge> edges = new HashSet<Edge>();

		BufferedReader reader = new BufferedReader(new FileReader(edgesFile));
		String line = reader.readLine();
		if (skipHeader)
		{
			line = reader.readLine();
		}
		while (line != null)
		{
			Edge e = EdgeIO.readFromTab(line.trim(), nodes);
			edges.add(e);
			line = reader.readLine();
		}

		reader.close();

		return edges;
	}

	/**
	 * Read all nodes from a file containing one node name per line
	 * 
	 * @param nodesFile the file containing the node data
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @param nodeAttributes the node attribute names
	 * 
	 * @return the set of nodes read from the file, or an empty set if no nodes were found
	 * @throws IOException when an error occurs during parsing
	 */
	public static Set<Node> readNodesFromFile(File nodesFile, NodeMapper nm, boolean skipHeader, List<String> nodeAttributes) throws IOException
	{
		Set<Node> nodes = new HashSet<Node>();

		BufferedReader reader = new BufferedReader(new FileReader(nodesFile));
		String line = reader.readLine();
		if (skipHeader)
		{
			line = reader.readLine();
		}
		while (line != null)
		{
			Node n = NodeIO.readFromTab(line.trim(), nodeAttributes);
			if (!nm.isContained(n, nodes))
			{
				nodes.add(n);
			}

			line = reader.readLine();
		}

		reader.close();
		return nodes;
	}

	/**
	 * Read the ID of a network from file. Specifically, a line of form "ID \t XYZ" is searched, and XYZ returned as the ID in integer form.
	 * In case more than one such line matches in the file, the first one is picked.
	 * 
	 * @param definitionFile the file containing the network definition data
	 * @return the ID of the network, as read from the file, or -1 if no ID was found or it could not be parsed as Integer
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static int readIDFromFile(File definitionFile) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(definitionFile));
		String line = reader.readLine();

		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			if (stok.nextToken().equals(ID_field))
			{
				String ID = stok.nextToken();
				reader.close();
				return Integer.parseInt(ID);
			}
			line = reader.readLine();
		}
		reader.close();
		return -1;
	}

	/**
	 * Read the name of a network from file. Specifically, a line of form "Name \t XYZ" is searched, and XYZ returned as the name.
	 * In case more than one such line matches in the file, the first one is picked.
	 * 
	 * @param definitionFile the file containing the network definition data
	 * @return the name of the network, as read from the file, or null if no name was found
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static String readNameFromFile(File definitionFile) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(definitionFile));
		String line = reader.readLine();

		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			if (stok.nextToken().equals(name_field))
			{
				String name = stok.nextToken();
				reader.close();
				return name;
			}
			line = reader.readLine();
		}
		reader.close();
		return null;
	}

	/**
	 * Read the type of a network from file. Specifically, a line of form "Type \t XYZ" is searched, and XYZ returned as the type.
	 * In case more than one such line matches in the file, the first one is picked.
	 * 
	 * @param definitionFile the file containing the network definition data
	 * @return the type of the network, as read from the file, or null if no type was found
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static String readTypeFromFile(File definitionFile) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(definitionFile));
		String line = reader.readLine();

		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			if (stok.nextToken().equals(type_field))
			{
				String type = stok.nextToken();
				reader.close();
				return type;
			}
			line = reader.readLine();
		}
		reader.close();
		return null;
	}

	/**
	 * Read the node attributes of a network from file. Specifically, a line of form "Attributes \t XYZ" is searched, and XYZ returned as the type.
	 * In case more than one such line matches in the file, the first one is picked.
	 * 
	 * @param definitionFile the file containing the network definition data
	 * @return the node attributes in the network, as read from the file, or an empty set if none were found
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static List<String> readAttributesFromFile(File definitionFile) throws IOException
	{
		List<String> attributes = new ArrayList<String>();

		BufferedReader reader = new BufferedReader(new FileReader(definitionFile));
		String line = reader.readLine();

		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			if (stok.nextToken().equals(attributes_field))
			{
				String attributeString = "";

				if (stok.hasMoreTokens())
				{
					attributeString = stok.nextToken();
				}

				StringTokenizer stok2 = new StringTokenizer(attributeString, ";");
				while (stok2.hasMoreTokens())
				{
					attributes.add(stok2.nextToken());
				}

				reader.close();
				return attributes;
			}
			line = reader.readLine();
		}
		reader.close();
		return attributes;
	}

}
