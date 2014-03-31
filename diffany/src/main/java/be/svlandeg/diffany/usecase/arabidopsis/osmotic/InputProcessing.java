package be.svlandeg.diffany.usecase.arabidopsis.osmotic;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import be.svlandeg.diffany.r.ExecuteR;

/**
 * This class reads and normalizes the raw input data with R scripts.
 * 
 * @author Sofie Van Landeghem
 */
public class InputProcessing
{

	/**
	 * Process the raw osmotic expression data with R.
	 * 
	 * Currently, the needed R script is loaded from the context, and is defined in the 'resources' folder of the Maven project.
	 * TODO v2.1: will this code work when packaged inside a jar or will we need to create a tmp file?
	 * 
	 * @throws URISyntaxException
	 * @throws IOException
	 */
	public void processOsmoticData(ExecuteR exeR, File osmoticStressDir) throws URISyntaxException, IOException
	{
		String path = osmoticStressDir.getAbsolutePath();
		System.out.println(" Reading " + path + ":");

		for (File f : osmoticStressDir.listFiles())
		{
			String fileName = f.getName();
			if (fileName.endsWith(".CEL"))
			{
				System.out.print("  " + fileName);
			}
		}
		System.out.println("");
		System.out.println("");

		String old_dir_path = exeR.changeExecutionDir(path);
		System.out.println(" Old WD in R: " + old_dir_path);
		System.out.println(" Set new WD in R to " + path);
		System.out.println("");

		// TODO V2.1: currently this assumes libs "affy", "affyPLM" and "org.Dm.eg.db" are pre-installed!
		URL script1URL = Thread.currentThread().getContextClassLoader().getResource("Rcode/ReadAffyData.R");
		System.out.println(" Executing script to read the raw data: " + script1URL);
		System.out.println(" (this may take a minute ... please be patient and do not interrupt the execution) ");
		exeR.executeScript(script1URL);
		System.out.println("");

		URL script2URL = Thread.currentThread().getContextClassLoader().getResource("Rcode/NormalizeAffyData.R");
		System.out.println(" Executing script to normalize the expression data: " + script2URL);
		System.out.println(" (this may take a minute ... please be patient and do not interrupt the execution) ");
		exeR.executeScript(script2URL);
		System.out.println("");
	}
}
