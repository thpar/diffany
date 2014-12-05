package be.svlandeg.diffany.core.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
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
	 * @param edgesFile the output file in which the edges will be written
	 * @param nodesFile the output file in which the nodes will be written
	 * @param definitionFile the output file in which the network definition will be written
	 * @param writeHeaders whether or not to write the headers of the nodes and edges files
	 * @see EdgeIO#writeToTab
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	protected static void writeNetworkToFiles(Network network, File edgesFile, File nodesFile, File definitionFile, boolean writeHeaders) throws IOException
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

		for (Edge e : network.getEdges())
		{
			edgeWriter.append(EdgeIO.writeToTab(e));
			edgeWriter.newLine();
			edgeWriter.flush();
		}

		edgeWriter.flush();
		edgeWriter.close();

		// NODES
		nodesFile.getParentFile().mkdirs();
		BufferedWriter nodeWriter = new BufferedWriter(new FileWriter(nodesFile));

		SortedSet<String> nodeAttributes = new TreeSet<String>();
		nodeAttributes.addAll(network.getAllNodeAttributes());

		if (writeHeaders)
		{
			nodeWriter.append(NodeIO.getHeader(nodeAttributes));
			nodeWriter.newLine();
			nodeWriter.flush();
		}

		// we assume the network has one node per unique ID!
		SortedMap<String, Node> sortedNodes = new TreeMap<String, Node>();
		for (Node n : network.getNodes())
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
		// TODO v2.2: create DefinitionIO to read and write the network definition
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
	 * @param dir the output dir in which the tab files will be written
	 * @param writeHeaders whether or not to write the headers of the nodes and edges files
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	public static void writeNetworkToDir(Network network, File dir, boolean writeHeaders) throws IOException
	{
		File edgeFile = new File(dir.getAbsolutePath() + "/" + default_edge_file);
		File nodeFile = new File(dir.getAbsolutePath() + "/" + default_node_file);
		File definitionFile = new File(dir.getAbsolutePath() + "/" + default_definition_file);
		writeNetworkToFiles(network, edgeFile, nodeFile, definitionFile, writeHeaders);

		if (network instanceof ConditionNetwork)
		{
			File default_conditions_File = new File(dir.getAbsolutePath() + "/" + default_conditions_file);
			writeConditionsToFiles((ConditionNetwork) network, default_conditions_File);
		}
	}

	/**
	 * Read a network from a directory in a Resource: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the dir in which the tab files are stored
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a Network representation of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	private static Network readInputNetworkFromResource(String dir, boolean skipHeader) throws IOException
	{
		InputStream edgeStream = NetworkIO.class.getResourceAsStream(dir + "/" + default_edge_file);
		InputStream nodeStream = NetworkIO.class.getResourceAsStream(dir + "/" + default_node_file);
		InputStream definitionStream = NetworkIO.class.getResourceAsStream(dir + "/" + default_definition_file);
		InputStream conditionsStream = NetworkIO.class.getResourceAsStream(dir + "/" + default_conditions_file);
		Network inputNetwork = readInputNetworkFromStreams(edgeStream, nodeStream, definitionStream, conditionsStream, skipHeader);
		edgeStream.close();
		nodeStream.close();
		definitionStream.close();
		if (conditionsStream != null)
		{
			conditionsStream.close();
		}
		return inputNetwork;
	}

	/**
	 * Read a network from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a Network representation of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	private static Network readInputNetworkFromDir(File dir, boolean skipHeader) throws IOException
	{
		InputStream edgeStream = new FileInputStream(new File(dir.getAbsolutePath() + "/" + default_edge_file));
		InputStream nodeStream = new FileInputStream(new File(dir.getAbsolutePath() + "/" + default_node_file));
		InputStream definitionStream = new FileInputStream(new File(dir.getAbsolutePath() + "/" + default_definition_file));
		InputStream conditionsStream = null;
		Network inputNetwork = null;
		try
		{
			conditionsStream = new FileInputStream(new File(dir.getAbsolutePath() + "/" + default_conditions_file));
		}
		catch (Exception e)
		{
			// it might be that there is no condition file
		}
		inputNetwork = readInputNetworkFromStreams(edgeStream, nodeStream, definitionStream, conditionsStream, skipHeader);
		edgeStream.close();
		nodeStream.close();
		definitionStream.close();
		if (conditionsStream != null)
		{
			conditionsStream.close();
		}

		return inputNetwork;
	}

	private static Network readInputNetworkFromStreams(InputStream edgeStream, InputStream nodeStream, InputStream definitionStream, InputStream conditionsStream, boolean skipHeader) throws IOException
	{

		Map<String, String> definitionMap = cacheDefinitionStream(definitionStream);
		int ID = readIDFromMap(definitionMap);
		String name = readNameFromMap(definitionMap);
		String type = readTypeFromMap(definitionMap);
		List<String> listedAttributes = readAttributesFromMap(definitionMap);
		Set<String> attributes = new HashSet<String>(listedAttributes);

		Set<Node> nodes = readNodesFromStream(nodeStream, skipHeader, listedAttributes);
		Set<Edge> edges = readEdgesFromStream(edgeStream, getMappedNodes(nodes), skipHeader);

		if (type.equals("ReferenceNetwork"))
		{
			ReferenceNetwork r = new ReferenceNetwork(name, ID, attributes);
			r.setNodesAndEdges(nodes, edges);
			return r;
		}

		else if (type.equals("ConditionNetwork"))
		{
			Set<Condition> conditions = readConditionsFromStream(conditionsStream);
			ConditionNetwork c = new ConditionNetwork(name, ID, attributes, conditions);
			c.setNodesAndEdges(nodes, edges);
			return c;
		}
		else if (type.equals("InputNetwork"))
		{
			InputNetwork c = new InputNetwork(name, ID, attributes);
			c.setNodesAndEdges(nodes, edges);
			return c;
		}

		throw new UnsupportedDataTypeException("Encountered unknown input network type: " + type);
	}

	/**
	 * Read a {@link ReferenceNetwork} from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a ReferenceNetwork representation of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static ReferenceNetwork readReferenceNetworkFromResource(String dir, boolean skipHeader) throws IOException
	{
		return (ReferenceNetwork) readInputNetworkFromResource(dir, skipHeader);
	}

	/**
	 * Read a {@link ReferenceNetwork} from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a ReferenceNetwork representation of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static ReferenceNetwork readReferenceNetworkFromDir(File dir, boolean skipHeader) throws IOException
	{
		return (ReferenceNetwork) readInputNetworkFromDir(dir, skipHeader);
	}

	/**
	 * Read a {@link ConditionNetwork} from a resource: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a ConditionNetwork representation of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static ConditionNetwork readConditionNetworkFromResource(String dir, boolean skipHeader) throws IOException
	{
		return (ConditionNetwork) readInputNetworkFromResource(dir, skipHeader);
	}

	/**
	 * Read a {@link ConditionNetwork} from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a ConditionNetwork representation of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static ConditionNetwork readConditionNetworkFromDir(File dir, boolean skipHeader) throws IOException
	{
		return (ConditionNetwork) readInputNetworkFromDir(dir, skipHeader);
	}

	/**
	 * Read an {@link InputNetwork} from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return an InputNetwork representation of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static InputNetwork readGenericInputNetworkFromDir(File dir, boolean skipHeader) throws IOException
	{
		return (InputNetwork) readInputNetworkFromDir(dir, skipHeader);
	}

	/**
	 * Read a set of {@link InputNetwork} from a directory, one network per subdirectory.
	 * 
	 * @param dir the output dir in which the subdirectories contain previously written tab files defining the different networks
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a set of InputNetwork representations of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static Set<InputNetwork> readGenericInputNetworksFromSubdirs(File dir, boolean skipHeader) throws IOException
	{
		Set<InputNetwork> networks = new HashSet<InputNetwork>();
		for (File f : dir.listFiles())
		{
			if (f.isDirectory())
			{
				InputNetwork net = readGenericInputNetworkFromDir(f, skipHeader);
				networks.add(net);
			}
		}
		return networks;
	}

	/**
	 * Read a network from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param reference the ReferenceNetwork linked to the read output network
	 * @param condNetworks the set of condition-specific networks linked to the read output network
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a Network representation of the nodes and edges in the files.
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	private static Network readOutputNetworkFromDir(File dir, ReferenceNetwork reference, Set<ConditionNetwork> condNetworks, boolean skipHeader) throws IOException
	{
		File edgeFile = new File(dir.getAbsolutePath() + "/" + default_edge_file);
		File nodeFile = new File(dir.getAbsolutePath() + "/" + default_node_file);
		File definitionFile = new File(dir.getAbsolutePath() + "/" + default_definition_file);

		int ID = readIDFromFile(definitionFile);
		String name = readNameFromFile(definitionFile);
		String type = readTypeFromFile(definitionFile);
		List<String> listedAttributes = readAttributesFromFile(definitionFile);

		Set<Node> nodes = readNodesFromFile(nodeFile, skipHeader, listedAttributes);
		Set<Edge> edges = readEdgesFromFile(edgeFile, getMappedNodes(nodes), skipHeader);

		if (type.equals("DifferentialNetwork"))
		{
			DifferentialNetwork d = new DifferentialNetwork(name, ID, reference, condNetworks);
			d.setNodesAndEdges(nodes, edges);
			return d;
		}
		else if (type.equals("ConsensusNetwork"))
		{
			Set<Network> allNetworks = new HashSet<Network>();
			allNetworks.add(reference);
			allNetworks.addAll(condNetworks);
			ConsensusNetwork o = new ConsensusNetwork(name, ID, allNetworks);
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
	 * @param reference the ReferenceNetwork linked to this differential network
	 * @param condNetworks the set of condition-specific networks linked to this differential network
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a DifferentialNetwork representation of the nodes and edges in the files
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static DifferentialNetwork readDifferentialNetworkFromDir(File dir, ReferenceNetwork reference, Set<ConditionNetwork> condNetworks, boolean skipHeader) throws IOException
	{
		return (DifferentialNetwork) readOutputNetworkFromDir(dir, reference, condNetworks, skipHeader);
	}

	/**
	 * Read a {@link ConsensusNetwork} from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param reference the ReferenceNetwork linked to this consensus network
	 * @param condNetworks the set of condition-specific networks linked to this consensus network
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return a ConsensusNetwork representation of the nodes and edges in the files.
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static ConsensusNetwork readConsensusNetworkFromDir(File dir, ReferenceNetwork reference, Set<ConditionNetwork> condNetworks, boolean skipHeader) throws IOException
	{
		return (ConsensusNetwork) readOutputNetworkFromDir(dir, reference, condNetworks, skipHeader);
	}

	/**
	 * Read all conditions from a stream containing one tab-delimited condition per line.
	 * 
	 * @param conditionsStream the stream containing the condition data
	 * @return the set of conditions read from the file, or an empty set if no conditions were found
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static Set<Condition> readConditionsFromStream(InputStream conditionsStream) throws IOException
	{
		Set<Condition> conditions = new HashSet<Condition>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(conditionsStream));
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

		return conditions;
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
		FileInputStream condStream = new FileInputStream(conditionsFile);
		Set<Condition> conds = readConditionsFromStream(condStream);
		condStream.close();
		return conds;
	}

	/**
	 * Read all edges from a stream containing one tab-delimited edge per line. As input, a set of nodes should be given, mapped by their unique IDs.
	 * 
	 * @param edgesStream the stream containing the edge data
	 * @param nodes the nodes relevant to the edges that will be read
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @return the set of edges read from the file, or an empty set if no edges were found
	 * @see EdgeIO#readFromTab
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static Set<Edge> readEdgesFromStream(InputStream edgesStream, Map<String, Node> nodes, boolean skipHeader) throws IOException
	{
		Set<Edge> edges = new HashSet<Edge>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(edgesStream));
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

		return edges;
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
		FileInputStream edgesStream = new FileInputStream(edgesFile);
		Set<Edge> edges = readEdgesFromStream(edgesStream, nodes, skipHeader);
		edgesStream.close();
		return edges;
	}

	/**
	 * Read all nodes from a stream containing one node name per line
	 * 
	 * @param nodesStream the stream containing the node data
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @param nodeAttributes the node attribute names
	 * 
	 * @return the set of nodes read from the file, or an empty set if no nodes were found
	 * @throws IOException when an error occurs during parsing
	 */
	public static Set<Node> readNodesFromStream(InputStream nodesStream, boolean skipHeader, List<String> nodeAttributes) throws IOException
	{
		Set<Node> nodes = new HashSet<Node>();
		Set<String> read_IDs = new HashSet<String>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(nodesStream));
		String line = reader.readLine();
		if (skipHeader)
		{
			line = reader.readLine();
		}
		while (line != null)
		{
			Node n = NodeIO.readFromTab(line.trim(), nodeAttributes);
			String ID = n.getID();
			if (!read_IDs.contains(ID))
			{
				read_IDs.add(ID);
				nodes.add(n);
			}
			line = reader.readLine();
		}

		return nodes;
	}

	/**
	 * Read all nodes from a file containing one node name per line
	 * 
	 * @param nodesFile the file containing the node data
	 * @param skipHeader whether or not the nodes and edges file contain a header
	 * @param nodeAttributes the node attribute names
	 * 
	 * @return the set of nodes read from the file, or an empty set if no nodes were found
	 * @throws IOException when an error occurs during parsing
	 */
	public static Set<Node> readNodesFromFile(File nodesFile, boolean skipHeader, List<String> nodeAttributes) throws IOException
	{
		FileInputStream nodesStream = new FileInputStream(nodesFile);
		Set<Node> nodes = readNodesFromStream(nodesStream, skipHeader, nodeAttributes);
		nodesStream.close();
		return nodes;
	}

	/**
	 * Read and parse a definition file from an input stream and store the key value pairs in a Map.
	 */
	private static Map<String, String> cacheDefinitionStream(InputStream definitionStream) throws IOException
	{
		BufferedReader reader = new BufferedReader(new InputStreamReader(definitionStream));
		String line = reader.readLine();

		Map<String, String> definitionMap = new HashMap<String, String>();

		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			String key = stok.nextToken();
			String value = new String();
			if (stok.hasMoreTokens())
			{
				value = stok.nextToken();
			}
			else
			{
				value = "";
			}
			definitionMap.put(key, value);
			line = reader.readLine();
		}
		return definitionMap;
	}

	/**
	 * Return the ID of a network from stream.
	 * 
	 * @param definitionMap key value pairs read from a definition file
	 * @return the ID of the network, as read from the file, or -1 if no ID was found or it could not be parsed as Integer
	 */
	public static int readIDFromMap(Map<String, String> definitionMap)
	{
		if (definitionMap.containsKey(ID_field))
		{
			return Integer.valueOf(definitionMap.get(ID_field));
		}
		return -1;
	}

	/**
	 * Read the ID of a network from stream. Specifically, a line of form "ID \t XYZ" is searched, and XYZ returned as the ID in integer form.
	 * In case more than one such line matches in the file, the first one is picked.
	 * 
	 * @param definitionStream the file containing the network definition data
	 * @return the ID of the network, as read from the file, or -1 if no ID was found or it could not be parsed as Integer
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static int readIDFromStream(InputStream definitionStream) throws IOException
	{
		BufferedReader reader = new BufferedReader(new InputStreamReader(definitionStream));
		String line = reader.readLine();

		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			if (stok.nextToken().equals(ID_field))
			{
				String ID = stok.nextToken();
				return Integer.parseInt(ID);
			}
			line = reader.readLine();
		}
		return -1;
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
		FileInputStream defStream = new FileInputStream(definitionFile);
		int id = readIDFromStream(defStream);
		defStream.close();
		return id;
	}

	/**
	 * Return the name of a network from a stream.
	 * 
	 * @param definitionMap key value pairs read from a definition file
	 * 
	 * @return the name of the network, as read from the stream, or null if no name was found
	 */
	private static String readNameFromMap(Map<String, String> definitionMap)
	{
		if (definitionMap.containsKey(name_field))
		{
			return definitionMap.get(name_field);
		}
		return null;
	}

	/**
	 * Read the name of a network from a stream. Specifically, a line of form "Name \t XYZ" is searched, and XYZ returned as the name.
	 * In case more than one such line matches in the file, the first one is picked.
	 * 
	 * @param definitionStream the stream containing the network definition data
	 * @return the name of the network, as read from the stream, or null if no name was found
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static String readNameFromStream(InputStream definitionStream) throws IOException
	{
		BufferedReader reader = new BufferedReader(new InputStreamReader(definitionStream));
		String line = reader.readLine();

		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			if (stok.nextToken().equals(name_field))
			{
				String name = stok.nextToken();
				return name;
			}
			line = reader.readLine();
		}
		return null;
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
		InputStream definitionStream = new FileInputStream(definitionFile);
		String name = readNameFromStream(definitionStream);
		definitionStream.close();
		return name;
	}

	private static String readTypeFromMap(Map<String, String> definitionMap)
	{
		if (definitionMap.containsKey(type_field))
		{
			return definitionMap.get(type_field);
		}
		return null;
	}

	/**
	 * Read the type of a network from stream. Specifically, a line of form "Type \t XYZ" is searched, and XYZ returned as the type.
	 * In case more than one such line matches in the file, the first one is picked.
	 * 
	 * @param definitionStream the stream containing the network definition data
	 * @return the type of the network, as read from the file, or null if no type was found
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static String readTypeFromStream(InputStream definitionStream) throws IOException
	{
		BufferedReader reader = new BufferedReader(new InputStreamReader(definitionStream));
		String line = reader.readLine();

		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			if (stok.nextToken().equals(type_field))
			{
				String type = stok.nextToken();
				return type;
			}
			line = reader.readLine();
		}
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
		FileInputStream defStream = new FileInputStream(definitionFile);
		String type = readTypeFromStream(defStream);
		defStream.close();
		return type;

	}

	/**
	 * Read the node attributes of a network from map.
	 * 
	 * @param definitionMap the hashmap containing the network definition data
	 * @return the node attributes in the network, as read from the file, or an empty set if none were found
	 */
	private static List<String> readAttributesFromMap(Map<String, String> definitionMap)
	{
		List<String> attributes = new ArrayList<String>();
		if (definitionMap.containsKey(attributes_field))
		{
			String attributeString = definitionMap.get(attributes_field);
			StringTokenizer stok = new StringTokenizer(attributeString, ";");
			while (stok.hasMoreTokens())
			{
				attributes.add(stok.nextToken());
			}
		}
		return attributes;
	}

	/**
	 * Read the node attributes of a network from stream. Specifically, a line of form "Attributes \t XYZ" is searched, and XYZ returned as the type.
	 * In case more than one such line matches in the file, the first one is picked.
	 * 
	 * @param definitionStream the stream containing the network definition data
	 * @return the node attributes in the network, as read from the file, or an empty set if none were found
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static List<String> readAttributesFromStream(InputStream definitionStream) throws IOException
	{
		List<String> attributes = new ArrayList<String>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(definitionStream));
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

				return attributes;
			}
			line = reader.readLine();
		}
		return attributes;
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
		FileInputStream defStream = new FileInputStream(definitionFile);
		List<String> attrs = readAttributesFromStream(defStream);
		defStream.close();
		return attrs;
	}
}
