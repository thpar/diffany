package be.svlandeg.diffany.cytoscape.layout;

import java.util.Set;

import org.cytoscape.model.CyNode;
import org.cytoscape.view.layout.AbstractLayoutAlgorithm;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.View;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.undo.UndoSupport;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.tasks.CopyLayoutTask;

/**
 * Layout algorithm that tries to layout a view identical to a given reference view, based on 
 * shared node names. 
 * 
 * @author Thomas Van Parys
 *
 */
public class CopyLayout extends AbstractLayoutAlgorithm{

	private Model model;
	private UndoSupport undo;
	private final static String NAME = "copylayout";

	public CopyLayout(Model model){
		super(NAME, "Copy Layout", null);
		UndoSupport undo = null;
		this.model = model;
		this.undo = undo;
	}
	

	@Override
	public TaskIterator createTaskIterator(CyNetworkView networkView,
			Object layoutContext, Set<View<CyNode>> nodesToLayOut,
			String layoutAttribute) {
		
		CopyLayoutTask task = new CopyLayoutTask(model.getCurrentProject(), 
				NAME, networkView, nodesToLayOut, layoutAttribute, undo);
		TaskIterator it = new TaskIterator();
		it.append(task);
		return it;
	}

}
