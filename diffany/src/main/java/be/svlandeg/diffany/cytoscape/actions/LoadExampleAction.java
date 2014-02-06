package be.svlandeg.diffany.cytoscape.actions;

import java.awt.event.ActionEvent;

import org.cytoscape.application.swing.AbstractCyAction;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.swing.DialogTaskManager;

import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.cytoscape.internal.Services;
import be.svlandeg.diffany.cytoscape.tasks.LoadExampleTaskFactory;

/**
 * Action that load {@link CyNetwork}s into a Cytoscape session, based on a given {@link Project}.
 * 
 * @author Thomas Van Parys
 *
 */
public class LoadExampleAction extends AbstractCyAction{

	private static final long serialVersionUID = 1L;
	private Services services;
	private Project exampleProject;
	private int runConfigurationID;
	
	/**
	 * Load an example {@link Project} into the current Cytoscape session. 
	 * Project settings are ignored. Only the input networks are read and translated to {@link CyNetwork}s.
	 * 
	 * @param services the Cytoscape {@link Services}s
	 * @param name display name in the menu
	 * @param exampleProject {@link Project} to load the source networks from
	 * @param runConfigurationID id of the example configuration to run in the project
	 */
	public LoadExampleAction(Services services, String name, Project exampleProject, int runConfigurationID) {
		super(name, services.getCyApplicationManager(), null, null);
		this.services = services;
		this.exampleProject = exampleProject;
		this.runConfigurationID = runConfigurationID;
		setPreferredMenu("Apps.Diffany.Examples");
	}


	@Override
	public void actionPerformed(ActionEvent e) {
		LoadExampleTaskFactory tf = new LoadExampleTaskFactory(services, exampleProject, runConfigurationID);

		if (tf.isReady()){
			TaskIterator it = tf.createTaskIterator();			
			DialogTaskManager dtm = services.getDialogTaskManager();
			dtm.execute(it);
		}

	}

}
