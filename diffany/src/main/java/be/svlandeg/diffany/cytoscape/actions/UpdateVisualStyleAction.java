package be.svlandeg.diffany.cytoscape.actions;

import java.awt.event.ActionEvent;

import org.cytoscape.application.swing.AbstractCyAction;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskManager;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.tasks.UpdateVisualStyleTask;
import be.svlandeg.diffany.cytoscape.tasks.UpdateVisualStyleTaskFactory;

/**
 * Run the {@link UpdateVisualStyleTask}
 * 
 * @author Thomas Van Parys
 *
 */
public class UpdateVisualStyleAction extends AbstractCyAction {
	
	private static final long serialVersionUID = -8382501915601785381L;
	private static final String BUTTON_TITLE = "Update Visual Style";
	private Model model;
	private CyProject cyProject;
	
	public UpdateVisualStyleAction(Model model, CyProject cyProject){
		super(BUTTON_TITLE);
		this.model = model;
		this.cyProject = cyProject;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		UpdateVisualStyleTaskFactory tf = new UpdateVisualStyleTaskFactory(model, cyProject);

		if (tf.isReady()){
			TaskIterator it = tf.createTaskIterator();			
			TaskManager tm = model.getServices().getTaskManager();
			tm.execute(it);
		}

	}

}
