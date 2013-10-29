package be.svlandeg.cytoscape.diffany.internal;

import java.util.Properties;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.service.util.AbstractCyActivator;
import org.osgi.framework.BundleContext;

/**
 * Defines a MenuAction and registers it as an OSGi service to the Cytoscape
 * application manager.
 * 
 * Generated from the cyaction-app Archetype.
 */
public class CyActivator extends AbstractCyActivator
{

	@Override
	public void start(BundleContext context) throws Exception
	{

		CyApplicationManager cyApplicationManager = getService(context, CyApplicationManager.class);

		MenuAction action = new MenuAction(cyApplicationManager, "Diffany");

		Properties properties = new Properties();

		registerAllServices(context, action, properties);
	}

}
