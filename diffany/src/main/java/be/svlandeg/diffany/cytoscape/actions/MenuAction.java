package be.svlandeg.diffany.cytoscape.actions;

import java.awt.event.ActionEvent;

import javax.swing.JOptionPane;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.application.swing.AbstractCyAction;

/**
 * Creates a new menu item in the Apps menu section.
 * 
 * Generated from the cyaction-app Archetype.
 * 
 * TODO: "Instead of using this class directly you should (strongly) consider
 * implementing a TaskFactory/Task pair. Doing so will allow your action to be
 * used outside of a Swing specific application (which the CyAction interface
 * binds you to)! "
 */
public class MenuAction extends AbstractCyAction
{

	/**
	 * Create this menu into the "Apps" folder of Cytoscape 3.X.X
	 */
	public MenuAction(CyApplicationManager cyApplicationManager, final String menuTitle)
	{

		super(menuTitle, cyApplicationManager, null, null);
		setPreferredMenu("Apps.Diffany");
	}

	/**
	 * Action called when the menu item of this app is selected.
	 * 
	 * TODO: create proper dialogue window, offering functionality to create,
	 * visualise and analyse 2 (or more) static networks + their differential
	 * network(s)
	 */
	public void actionPerformed(ActionEvent e)
	{
		// Write your own function here.
		// git test
		// git test 2
		// git test 3
		JOptionPane.showMessageDialog(null, "Something differentially :-)");
	}
}
