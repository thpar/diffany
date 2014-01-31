package be.svlandeg.diffany.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeSet;

import javax.activation.UnsupportedDataTypeException;

import be.svlandeg.diffany.concepts.Condition;
import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Node;
import be.svlandeg.diffany.concepts.OverlappingNetwork;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * This class allows reading or writing a {@link Network} from File.
 * 
 * @author Sofie Van Landeghem
 */
public class NetworkIO
{

	private static String default_conditions_file = "conditions.tab";
	private static String default_edge_file = "edges.tab";
	private static String default_node_file = "nodes.tab";
	private static String default_definition_file = "network.tab";
	
	private static String name_field = "Name";
	private static String type_field = "Type";
	

	/**
	 * Write the string representation of the network: All edges to one File, and all nodes to another, as specified by the parameters.
	 * Additionally, the network name and type are written to a defitionsfile
	 * 
	 * @param network the {@link Network} that needs to be written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param edgesFile the output file in which the edges will be written
	 * @param nodesFile the output file in which the nodes will be written
	 * @param definitionFile the output file in which the network definition will be written
	 * @see EdgeIO#writeToTab
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	protected static void writeNetworkToFiles(Network network, NodeMapper nm, File edgesFile, File nodesFile, File definitionFile) throws IOException
	{
		// EDGES
		edgesFile.getParentFile().mkdirs();
		BufferedWriter edgeWriter = new BufferedWriter(new FileWriter(edgesFile));
		for (Edge e : network.getEdges())
		{
			edgeWriter.append(EdgeIO.writeToTab(e));
			edgeWriter.newLine();
			edgeWriter.flush();
		}
		edgeWriter.flush();
		edgeWriter.close();

		// NODES
		// TODO v2.0: this code now simply prints the node name, but could be made more general through NodeIO and Node IDs?
		nodesFile.getParentFile().mkdirs();
		BufferedWriter nodeWriter = new BufferedWriter(new FileWriter(nodesFile));
		
		Set<Node> unduplicatedNodes = new HashSet<Node>();
		for (Node n : network.getNodes())
		{
			if (! nm.isContained(n, unduplicatedNodes))
			{
				unduplicatedNodes.add(n);
			}
		}
		
		TreeSet<String> sortedNodes = new TreeSet<String>();
		for (Node n : unduplicatedNodes)
		{
			sortedNodes.add(n.getName());
		}
		for (String name : sortedNodes)
		{
			nodeWriter.append(name);
			nodeWriter.newLine();
			nodeWriter.flush();
		}
		nodeWriter.flush();
		nodeWriter.close();
		
		// DEFINITION: NAME and TYPE (CLASS)
		definitionFile.getParentFile().mkdirs();
		BufferedWriter defWriter = new BufferedWriter(new FileWriter(definitionFile));
		
		defWriter.append(name_field + "\t" + network.getName());
		defWriter.newLine();
		
		String networkClass = network.getClass().getSimpleName();
		defWriter.append(type_field + "\t" + networkClass);
		defWriter.newLine();
		
		defWriter.flush();
		defWriter.close();
	}
	
	/**
	 * Write the string representation of the reference network to (subfiles of) a directory.
	 * 
	 * @param network the {@link ReferenceNetwork} that needs to be written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param dir the output dir in which the tab files will be written
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	public static void writeReferenceNetworkToDir(ReferenceNetwork network, NodeMapper nm, File dir) throws IOException
	{
		File edgeFile = new File(dir.getAbsolutePath() + "/" + default_edge_file);
		File nodeFile = new File(dir.getAbsolutePath() + "/" + default_node_file);
		File definitionFile = new File(dir.getAbsolutePath() + "/" + default_definition_file);
		writeReferenceNetworkToFiles(network, nm, edgeFile, nodeFile, definitionFile);
	}
	
	/**
	 * Write the string representation of the reference network: 
	 * Different files will contain the edges information, nodes and network definition.
	 * 
	 * @param network the {@link ReferenceNetwork} that needs to be written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param edgesFile the output file in which the edges will be written
	 * @param nodesFile the output file in which the nodes will be written
	 * @param definitionFile the output file in which the network definition will be written
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	public static void writeReferenceNetworkToFiles(ReferenceNetwork network, NodeMapper nm, File edgesFile, File nodesFile, File definitionFile)
			throws IOException
	{
		writeNetworkToFiles(network, nm, edgesFile, nodesFile, definitionFile);
	}
	
	
	/**
	 * Write the string representation of the differential network to (subfiles of) a directory.
	 * 
	 * @param network the {@link DifferentialNetwork} that needs to be written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param dir the output dir in which the tab files will be written
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	public static void writeDifferentialNetworkToDir(DifferentialNetwork network, NodeMapper nm, File dir) throws IOException
	{
		File edgeFile = new File(dir.getAbsolutePath() + "/" + default_edge_file);
		File nodeFile = new File(dir.getAbsolutePath() + "/" + default_node_file);
		File definitionFile = new File(dir.getAbsolutePath() + "/" + default_definition_file);
		writeDifferentialNetworkToFiles(network, nm, edgeFile, nodeFile, definitionFile);
	}
	
	/**
	 * Write the string representation of the differential network: 
	 * Different files will contain the edges information, nodes and network definition.
	 * 
	 * @param network the {@link DifferentialNetwork} that needs to be written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param edgesFile the output file in which the edges will be written
	 * @param nodesFile the output file in which the nodes will be written
	 * @param definitionFile the output file in which the network definition will be written
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	public static void writeDifferentialNetworkToFiles(DifferentialNetwork network, NodeMapper nm, File edgesFile, File nodesFile, File definitionFile)
			throws IOException
	{
		writeNetworkToFiles(network, nm, edgesFile, nodesFile, definitionFile);
	}
	
	
	
	/**
	 * Write the string representation of the overlapping network to (subfiles of) a directory.
	 * 
	 * @param network the {@link OverlappingNetwork} that needs to be written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param dir the output dir in which the tab files will be written
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	public static void writeOverlappingNetworkToDir(OverlappingNetwork network, NodeMapper nm, File dir) throws IOException
	{
		File edgeFile = new File(dir.getAbsolutePath() + "/" + default_edge_file);
		File nodeFile = new File(dir.getAbsolutePath() + "/" + default_node_file);
		File definitionFile = new File(dir.getAbsolutePath() + "/" + default_definition_file);
		writeOverlappingNetworkToFiles(network, nm, edgeFile, nodeFile, definitionFile);
	}
	
	/**
	 * Write the string representation of the overlapping network: 
	 * Different files will contain the edges information, nodes and network definition.
	 * 
	 * @param network the {@link OverlappingNetwork} that needs to be written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param edgesFile the output file in which the edges will be written
	 * @param nodesFile the output file in which the nodes will be written
	 * @param definitionFile the output file in which the network definition will be written
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	public static void writeOverlappingNetworkToFiles(OverlappingNetwork network, NodeMapper nm, File edgesFile, File nodesFile, File definitionFile)
			throws IOException
	{
		writeNetworkToFiles(network, nm ,edgesFile, nodesFile, definitionFile);
	}
	
	
	/**
	 * Write the string representation of the condition-specific network to (subfiles of) a directory.
	 * 
	 * @param network the {@link ConditionNetwork} that needs to be written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param dir the output dir in which the tab files will be written
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	public static void writeConditionNetworkToDir(ConditionNetwork network, NodeMapper nm, File dir) throws IOException
	{
		File edgeFile = new File(dir.getAbsolutePath() + "/" + default_edge_file);
		File nodeFile = new File(dir.getAbsolutePath() + "/" + default_node_file);
		File definitionFile = new File(dir.getAbsolutePath() + "/" + default_definition_file);
		File default_conditions_File = new File(dir.getAbsolutePath() + "/" + default_conditions_file);
		writeConditionNetworkToFiles(network, nm, edgeFile, nodeFile, definitionFile, default_conditions_File);
	}
	
	/**
	 * Write the string representation of the condition-specific network: 
	 * Different files will contain the edges information, nodes, network definition and conditions.
	 * 
	 * @param network the {@link ConditionNetwork} that needs to be written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param edgesFile the output file in which the edges will be written
	 * @param nodesFile the output file in which the nodes will be written
	 * @param definitionFile the output file in which the network definition will be written
	 * @param network the Network that needs to be written
	 * 
	 * @throws IOException when an error occurs during writing
	 */
	public static void writeConditionNetworkToFiles(ConditionNetwork network, NodeMapper nm, File edgesFile, File nodesFile, File definitionFile, File conditionsFile)
			throws IOException
	{
		writeNetworkToFiles(network, nm, edgesFile, nodesFile, definitionFile);
		
		Set<Condition> conditions = network.getConditions();
	
		conditionsFile.getParentFile().mkdirs();
		BufferedWriter writer = new BufferedWriter(new FileWriter(conditionsFile));
		for (Condition c : conditions)
		{
			writer.append(c.getDescription());
			for (String ont: c.getOntologyTerms())
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
	 * Read a network from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @return a Network representation of the nodes and edges in the files.
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	private static Network readInputNetworkFromDir(File dir, NodeMapper nm) throws IOException
	{
		File edgeFile = new File(dir.getAbsolutePath() + "/" + default_edge_file);
		File nodeFile = new File(dir.getAbsolutePath() + "/" + default_node_file);
		File definitionFile = new File(dir.getAbsolutePath() + "/" + default_definition_file);
		File conditionsFile = new File(dir.getAbsolutePath() + "/" + default_conditions_file);
		
		Set<Edge> edges = readEdgesFromFile(edgeFile);
		Set<Node> nodes = readNodesFromFile(nodeFile, nm);
		
		String name = readNameFromFile(definitionFile);
		String type = readTypeFromFile(definitionFile);
		
		if (type.equals("ReferenceNetwork"))
		{
			ReferenceNetwork r = new ReferenceNetwork(name, nm);
			r.setNodesAndEdges(nodes, edges);
			return r;
		}
		
		else if (type.equals("ConditionNetwork"))
		{
			Set<Condition> conditions = readConditionsFromFile(conditionsFile);
			ConditionNetwork c = new ConditionNetwork(name, conditions, nm);
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
	 * @return a ReferenceNetwork representation of the nodes and edges in the files.
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static ReferenceNetwork readReferenceNetworkFromDir(File dir, NodeMapper nm) throws IOException
	{
		return (ReferenceNetwork) readInputNetworkFromDir(dir, nm);
	}
	
	/**
	 * Read a {@link ConditionNetwork} from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @return a ConditionNetwork representation of the nodes and edges in the files.
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static ConditionNetwork readConditionNetworkFromDir(File dir, NodeMapper nm) throws IOException
	{
		return (ConditionNetwork) readInputNetworkFromDir(dir, nm);
	}
	
	/**
	 * Read a network from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param reference the ReferenceNetwork linked to the read output network
	 * @param condNetworks the set of condition-specific networks linked to the read output network
	 * @return a Network representation of the nodes and edges in the files.
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	private static Network readOutputNetworkFromDir(File dir, NodeMapper nm, ReferenceNetwork reference, Set<ConditionNetwork> condNetworks) throws IOException
	{
		File edgeFile = new File(dir.getAbsolutePath() + "/" + default_edge_file);
		File nodeFile = new File(dir.getAbsolutePath() + "/" + default_node_file);
		File definitionFile = new File(dir.getAbsolutePath() + "/" + default_definition_file);
		
		Set<Edge> edges = readEdgesFromFile(edgeFile);
		Set<Node> nodes = readNodesFromFile(nodeFile, nm);
		
		String name = readNameFromFile(definitionFile);
		String type = readTypeFromFile(definitionFile);
		
		if (type.equals("DifferentialNetwork"))
		{
			DifferentialNetwork d = new DifferentialNetwork(name, reference, condNetworks, nm);
			d.setNodesAndEdges(nodes, edges);
			return d;
		}
		else if (type.equals("OverlappingNetwork"))
		{
			Set<Network> allNetworks = new HashSet<Network>();
			allNetworks.add(reference);
			allNetworks.addAll(condNetworks);
			OverlappingNetwork o = new OverlappingNetwork(name, allNetworks, nm);
			o.setNodesAndEdges(nodes, edges);
			return o;
		}
		
		throw new UnsupportedDataTypeException("Encountered unknown output network type: " + type);
	}
	
	/**
	 * Read a {@link DifferentialNetwork} from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param reference the ReferenceNetwork linked to this differential network
	 * @param condNetworks the set of condition-specific networks linked to this differential network
	 * @return a DifferentialNetwork representation of the nodes and edges in the files
	 *
	 * @throws IOException when an error occurs during reading
	 */
	public static DifferentialNetwork readDifferentialNetworkFromDir(File dir, NodeMapper nm, ReferenceNetwork reference, Set<ConditionNetwork> condNetworks) throws IOException
	{
		return (DifferentialNetwork) readOutputNetworkFromDir(dir, nm, reference, condNetworks);
	}
	
	/**
	 * Read a {@link OverlappingNetwork} from a directory: all edges from one File (edges.tab), and all nodes from another (nodes.tab).
	 * 
	 * @param dir the output dir in which the tab files were previously written
	 * @param nm the {@link NodeMapper} object that determines equality between nodes
	 * @param reference the ReferenceNetwork linked to this overlapping network
	 * @param condNetworks the set of condition-specific networks linked to this overlapping network
	 * @return a OverlappingNetwork representation of the nodes and edges in the files.
	 * 
	 * @throws IOException when an error occurs during reading
	 */
	public static OverlappingNetwork readOverlappingNetworkFromDir(File dir, NodeMapper nm, ReferenceNetwork reference, Set<ConditionNetwork> condNetworks) throws IOException
	{
		return (OverlappingNetwork) readOutputNetworkFromDir(dir, nm, reference, condNetworks);
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
	 * Read all edges from a file containing one tab-delimited edge per line
	 * 
	 * @param edgesFile the file containing the edge data
	 * @return the set of edges read from the file, or an empty set if no edges were found
	 * @see EdgeIO#readFromTab
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static Set<Edge> readEdgesFromFile(File edgesFile) throws IOException
	{
		Set<Edge> edges = new HashSet<Edge>();

		BufferedReader reader = new BufferedReader(new FileReader(edgesFile));
		String line = reader.readLine();
		while (line != null)
		{
			Edge e = EdgeIO.readFromTab(line);
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
	 * @return the set of nodes read from the file, or an empty set if no nodes were found
	 * 
	 * @throws IOException when an error occurs during parsing
	 */
	public static Set<Node> readNodesFromFile(File nodesFile, NodeMapper nm) throws IOException
	{
		Set<Node> nodes = new HashSet<Node>();

		BufferedReader reader = new BufferedReader(new FileReader(nodesFile));
		String line = reader.readLine();
		
		while (line != null)
		{
			// TODO v2.0: this code now simply reads the node name, but could be made more general through NodeIO and Node IDs?
			Node n = new Node(line.trim());
			if (! nm.isContained(n, nodes))
			{
				nodes.add(n);
			}
			
			line = reader.readLine();
		}

		reader.close();
		return nodes;
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
			StringTokenizer stok = new StringTokenizer(line,"\t");
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
			StringTokenizer stok = new StringTokenizer(line,"\t");
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

}
