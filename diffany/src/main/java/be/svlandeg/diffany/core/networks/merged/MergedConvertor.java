package be.svlandeg.diffany.core.networks.merged;

import java.util.HashSet;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.networks.OverlappingNetwork;
import be.svlandeg.diffany.core.semantics.NodeMapper;

/**
 * This class allows quick conversion between merged input networks with condition edges, and normal input networks 
 * with normal edges.
 * 
 * @author Sofie Van Landeghem
 */
public class MergedConvertor
{

	/**
	 * Convert a set of condition edges to a set of normal edges.
	 * @param cEdges the condition-specific edges
	 * @return the set of normal edges
	 */
	public static Set<Edge> castToNormalEdges(Set<MergedEdge> cEdges)
	{
		Set<Edge> edges = new HashSet<Edge>();
		for (MergedEdge e : cEdges)
		{
			edges.add(e);
		}
		return edges;
	}

	/**
	 * Convert a set of normal edges to a set of condition edges.
	 * @param edges a set of normal edges which are actually condition-specific (and can thus be cast!)
	 * @return the set of condition-specific edges
	 * @throws ClassCastException when the set of input edges are not of the type MergedEdge
	 */
	public static Set<MergedEdge> castToConditionEdges(Set<Edge> edges)
	{
		Set<MergedEdge> cEdges = new HashSet<MergedEdge>();
		for (Edge e : edges)
		{
			cEdges.add((MergedEdge) e);
		}
		return cEdges;
	}
	
	/**
	 * Convert a set of input networks to one large merged network
	 * @param inputs the original set of separate input networks (assumed not-null and not empty!)
	 * @return the merged input network
	 */
	@SuppressWarnings("unused")
	public static MergedInputNetwork convertInput(Set<InputNetwork> inputs)
	{
		// TODO
		return null;
	}

	/**
	 * Convert a set of differential networks to one large merged network
	 * The name is created by appending all individual input names to the prefix 'merged_' and the ID is determined as the maximum + 1.
	 * 
	 * @param diffSet the original set of separate differential networks (assumed not-null and not empty!)
	 * @return the merged differential network
	 */
	public static MergedDifferentialNetwork convertDifferentials(Set<DifferentialNetwork> diffSet)
	{
		Set<InputNetwork> inputs = new HashSet<InputNetwork>();
		SortedSet<String> names = new TreeSet<String>();
		SortedSet<Integer> IDs = new TreeSet<Integer>();
		Set<Node> nodes = new HashSet<Node>();
		Set<MergedEdge> edges = new HashSet<MergedEdge>();
		NodeMapper nm = null;

		for (DifferentialNetwork dn : diffSet)
		{
			// all names are stored alphabetically, IDs as well
			names.add(dn.getName());
			IDs.add(dn.getID());
			
			// all input networks
			inputs.add(dn.getReferenceNetwork());
			inputs.addAll(dn.getConditionNetworks());

			// the node mapper object should be the same one across all separate networks
			if (nm == null)
			{
				nm = dn.getNodeMapper();
			}
			else if (!nm.equals(dn.getNodeMapper()))
			{
				String errormsg = "The differential output network(s) do not have the same nodemapper instance!";
				throw new RuntimeException(errormsg);
			}
			
			// add new nodes that have not been added to the set previously
			Set<Node> theseNodes = dn.getNodes();
			for (Node n : theseNodes)
			{
				if (! nm.isContained(n, nodes))
				{
					nodes.add(n);
				}
			}
			
			// TODO: edges conversion
		}

		String name = "merged_";
		for (String n : names)
		{
			name += n + "_";
		}
		
		int ID = IDs.last() + 1;

		MergedInputNetwork input = convertInput(inputs);
		MergedDifferentialNetwork result = new MergedDifferentialNetwork(name, ID, nodes, edges, nm, input);
		
		return result;
	}
	
	/**
	 * Convert a set of overlapping networks to one large merged network
	 * @param overlapSet the original set of separate overlapping networks (assumed not-null and not empty!)
	 * @return the merged overlapping network
	 */
	@SuppressWarnings("unused")
    public static MergedOverlappingNetwork convertOverlapping(Set<OverlappingNetwork> overlapSet)
	{
		// TODO
		return null;
	}

}
