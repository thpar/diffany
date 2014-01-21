package be.svlandeg.diffany.cytoscape.layout;

import java.util.Set;

import org.cytoscape.model.CyNode;
import org.cytoscape.view.layout.AbstractLayoutAlgorithm;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.View;
import org.cytoscape.work.TaskIterator;

/**
 * Layout algorithm that tries to layout a view identical to a given reference view, based on 
 * shared node names. 
 * 
 * @author Thomas Van Parys
 *
 */
public class CopyLayout extends AbstractLayoutAlgorithm{

	private CyNetworkView refView;


	public CopyLayout(CyNetworkView refView){
		super("copylayout", "Copy Layout", null);
		this.refView = refView;
	}
	

	@Override
	public TaskIterator createTaskIterator(CyNetworkView networkView,
			Object layoutContext, Set<View<CyNode>> nodesToLayOut,
			String layoutAttribute) {
		// TODO Auto-generated method stub
		return null;
	}

}
