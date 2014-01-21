package be.svlandeg.diffany.cytoscape.tasks;

import java.util.Set;

import org.cytoscape.model.CyNode;
import org.cytoscape.view.layout.AbstractLayoutTask;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.View;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.undo.UndoSupport;

public class CopyLayoutTask extends AbstractLayoutTask{

	public CopyLayoutTask(String displayName, CyNetworkView networkView,
			Set<View<CyNode>> nodesToLayOut, String layoutAttribute,
			UndoSupport undo) {
		super(displayName, networkView, nodesToLayOut, layoutAttribute, undo);
		// TODO Auto-generated constructor stub
	}

	@Override
	protected void doLayout(TaskMonitor taskMonitor) {
		// TODO Auto-generated method stub
		
	}

}
