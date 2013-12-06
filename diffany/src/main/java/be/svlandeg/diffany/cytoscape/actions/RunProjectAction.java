package be.svlandeg.diffany.cytoscape.actions;

import java.awt.event.ActionEvent;

import org.cytoscape.application.swing.AbstractCyAction;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.swing.DialogTaskManager;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.tasks.RunProjectTaskFactory;
import be.svlandeg.diffany.cytoscape.tasks.TestTaskFactory;

public class RunProjectAction extends AbstractCyAction {

	private Model model;

	public RunProjectAction(Model model, String menuTitle) {
		super(menuTitle, model.getServices().getCyApplicationManager(), null, null);
		setPreferredMenu("Apps.Diffany");
		this.model = model;
	}

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	@Override
	public void actionPerformed(ActionEvent e) {
		RunProjectTaskFactory tf = new RunProjectTaskFactory(model);
		
		if (tf.isReady()){
			TaskIterator it = tf.createTaskIterator();			
			DialogTaskManager dtm = model.getServices().getDialogTaskManager();
			dtm.execute(it);
		}
	}

}
