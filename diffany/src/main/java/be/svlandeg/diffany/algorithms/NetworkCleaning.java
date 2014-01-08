package be.svlandeg.diffany.algorithms;

import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.concepts.*;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * This class provides generic methods useful for network cleaning before or
 * after applying differential algorithms.
 * 
 * @author Sofie Van Landeghem
 */
public class NetworkCleaning
{
	

	/**
	 * Clean a network by adjusting symmetrical edges between two nodes that
	 * also have a directed node between them. In this case, the symmetrical
	 * edge is only kept in the direction that otherwise does not have an edge.
	 * 
	 * @param net the network that needs cleaning
	 **/
	public void directSymmetricalWhenOverlapping(DifferentialNetwork net, EdgeOntology eo)
	{
		/*
		Set<Edge> copySet = new HashSet<Edge>();	// needed to avoid ConcurrentModificationException
		copySet.addAll(net.getEdges());
		for (Edge symmE : copySet)
		{
			// retrieve for each symmetrical edge whether there are directed
			// overlapping edges
			if (symmE.isSymmetrical())
			{
				Node source = symmE.getSource();
				Node target = symmE.getTarget();
				String symmType = symmE.getType();

				Edge fwdDirection = new Edge(source, target, symmE);
				fwdDirection.makeSymmetrical(false);
				
				fwdDirection.setType(getDirectedType(eo, symmType));
				
				Set<Edge> others = net.getAllEdges(source, target);
				for (Edge fwdE : others)
				{
					if (!symmE.equals(fwdE) && !fwdE.isSymmetrical())
					{
						fwdDirection = null;
					}
				}

				Edge reverseDirection = new Edge(target, source, symmE);
				reverseDirection.makeSymmetrical(false);
				reverseDirection.setType(getDirectedType(eo, symmE.getType()));
				
				others = net.getAllEdges(target, source);
				for (Edge backE : others)
				{
					if (!symmE.equals(backE) && !backE.isSymmetrical())
					{
						reverseDirection = null;
					}
				}
				// keep only the reverse direction
				if (fwdDirection == null && reverseDirection != null)
				{
					net.removeEdge(symmE);
					net.addEdge(reverseDirection);
				}
				// keep only the forward direction
				if (reverseDirection == null && fwdDirection != null)
				{
					net.removeEdge(symmE);
					net.addEdge(fwdDirection);
				}
				// keep none of the directions
				if (reverseDirection == null && fwdDirection == null)
				{
					net.removeEdge(symmE);
				}
			}
		}*/
	}
	

	/**
	 * Remove edges in the network that are symmetrical and are represented
	 * twice (source-target and target-source). One of the two is removed only
	 * then when the type, weight and negation are all equal.
	 * 
	 * @param net the network that needs cleaning
	 * @param nm the node mapper for the network
	 **/
	public void removeRedundantEdges(Network net, NodeMapper nm)
	{
		Set<Edge> removed_edges = new HashSet<Edge>();
		
		// remove duplicate symmetrical edges between source-target and
		// target-source
		for (Node n1 : net.getNodes())
		{
			String name1 = n1.getName();
			for (Node n2 : net.getNodes())
			{
				String name2 = n2.getName();
				if (name1.compareTo(name2) < 0)
				{
					Set<Edge> all_edges = net.getAllEdges(n1, n2, nm);
					for (Edge et : all_edges)
					{
						for (Edge eb : all_edges)
						{
							if (!et.equals(eb) && !removed_edges.contains(et))
							{
								if ((et.isSymmetrical() && eb.isSymmetrical()) && (et.getType().equals(eb.getType())))
								{
									if ((et.getWeight() == eb.getWeight()) && (et.isNegated() == eb.isNegated()))
									{
										net.removeEdge(eb);
										removed_edges.add(eb);
									}
								}
							}
						}
					}
				}
			}
		}
	}

}
