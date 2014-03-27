package be.svlandeg.diffany.core.visualstyle;

import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;


/**
 * This abstract class tells the GUI how to draw certain edges, depending on definitions from an {@link TreeEdgeOntology}.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class TreeEdgeDrawing extends EdgeDrawing
{

	protected TreeEdgeOntology teo;
	
	/**
	 * Create a new EdgeDrawing object, given a specific TreeEdgeOntology object.
	 * @param teo the TreeEdgeOntology that can semantically interpret the interaction types through a tree structure
	 */
	public TreeEdgeDrawing(TreeEdgeOntology teo)
	{
		this.teo = teo;
	}

}
