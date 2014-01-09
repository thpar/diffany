package be.svlandeg.diffany.cytoscape.actions;

import java.awt.event.ActionEvent;

import org.cytoscape.application.swing.AbstractCyAction;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.swing.DialogTaskManager;

import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.cytoscape.tasks.LoadExampleTaskFactory;
import be.svlandeg.diffany.internal.Services;

public class LoadExampleAction extends AbstractCyAction{

	private Services services;
	private Project exampleProject;
	
	
	public LoadExampleAction(Services services, String name, Project exampleProject) {
		super(name, services.getCyApplicationManager(), null, null);
		this.services = services;
		this.exampleProject = exampleProject;
		setPreferredMenu("Apps.Diffany.Examples");
	}


	@Override
	public void actionPerformed(ActionEvent e) {
		LoadExampleTaskFactory tf = new LoadExampleTaskFactory(services, exampleProject);

		if (tf.isReady()){
			TaskIterator it = tf.createTaskIterator();			
			DialogTaskManager dtm = services.getDialogTaskManager();
			dtm.execute(it);
		}

	}

}
