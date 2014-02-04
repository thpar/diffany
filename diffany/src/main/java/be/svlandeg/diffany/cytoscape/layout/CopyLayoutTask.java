package be.svlandeg.diffany.cytoscape.layout;

import java.util.Set;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;
import org.cytoscape.view.layout.AbstractLayoutTask;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.View;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.undo.UndoSupport;

import be.svlandeg.diffany.cytoscape.CyProject;

/**
 * Layout task that uses a given {@link CyProject} to synchronize a given view with its reference network. 
 * Synchronization is attempted based on node names.
 * 
 * @author Thomas Van Parys
 *
 */
public class CopyLayoutTask extends AbstractLayoutTask{

	private CyProject cyProject;

	public CopyLayoutTask(CyProject cyProject, String displayName, CyNetworkView networkView,
			Set<View<CyNode>> nodesToLayOut, String layoutAttribute,
			UndoSupport undo) {
		super(displayName, networkView, nodesToLayOut, layoutAttribute, undo);
		this.cyProject = cyProject;
	}

	@Override
	protected void doLayout(TaskMonitor taskMonitor) {
		CyNetwork ref = cyProject.getReferenceNetwork();
		//TODO get views...
	}

}
