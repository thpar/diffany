package be.svlandeg.diffany.io;

import be.svlandeg.diffany.concepts.Logger;
import be.svlandeg.diffany.concepts.Project;

/**
 * This class allows reading or writing a {@link Project} from File.
 * 
 * @author Sofie Van Landeghem
 */
public class ProjectIO
{
	
	/**
	 * Save project data to a specific file location.
	 * 
	 * @param p the project that will be saved
	 * @param fileLocation the location where the project should be saved
	 */
	public static void saveProject(Project p, String fileLocation)
	{
		//TODO v1.1: implement!
		throw new UnsupportedOperationException("Saving of project not yet implemented");
	}
	

	/**
	 * Load project data from a specific file location. 
	 * Make sure all restrictions on number of networks are respected during the load!
	 *
	 * @param fileLocation the location from where the project should be loaded
	 */
	public static Project loadFromFile(String fileLocation)
	{
		//TODO v1.1: implement!
		Logger logger = new Logger();
		throw new UnsupportedOperationException("Loading of project not yet implemented");
	}

}
