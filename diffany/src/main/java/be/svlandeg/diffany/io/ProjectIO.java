package be.svlandeg.diffany.io;


import be.svlandeg.diffany.concepts.Project;

/**
 * This class allows reading or writing a {@link Project} from File.
 * Will be implemented in v.2.0, currently this functionality is not available.
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
	@SuppressWarnings("unused")
    private static void saveProject(Project p, String fileLocation)
	{
		//TODO v2.0: implement!
		throw new UnsupportedOperationException("Saving of project not yet implemented");
	}
	

	/**
	 * Load project data from a specific file location. 
	 * Make sure all restrictions on number of networks are respected during the load!
	 *
	 * @param fileLocation the location from where the project should be loaded
	 */
	@SuppressWarnings("unused")
    private static Project loadFromFile(String fileLocation)
	{
		//TODO v2.0: implement!
		throw new UnsupportedOperationException("Loading of project not yet implemented");
	}

}
