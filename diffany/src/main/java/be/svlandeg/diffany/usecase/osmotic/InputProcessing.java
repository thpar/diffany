package be.svlandeg.diffany.usecase.osmotic;

import java.io.File;
import java.io.IOException;
import java.net.URL;

import be.svlandeg.diffany.r.ExecuteR;

/**
 * This class reads and processes the raw input data.
 * 
 * @author Sofie Van Landeghem
 */
public class InputProcessing
{

	/**
	 * TODO documentation v2.0
	 * 
	 * Currently, the needed R script is loaded from the context, and is defined in the 'resources' folder
	 * of the Maven project.
	 */
	public void processOsmoticData(ExecuteR exeR, File osmoticStressDir) throws IOException
	{
		String path = osmoticStressDir.getAbsolutePath();
		System.out.println("reading " + path + ":");
		
		for (File f : osmoticStressDir.listFiles())
		{
			String fileName = f.getName();
			if (fileName.endsWith(".CEL"))
			{
				System.out.print(" " + fileName);
			}
		}
		System.out.println("");
		System.out.println("");
		
		String old_dir_path = exeR.changeExecutionDir(path);
		System.out.println("Old WD in R: " + old_dir_path);
		System.out.println("Set new WD in R to " + path);
		System.out.println("");
		
		URL scriptURL = Thread.currentThread().getContextClassLoader().getResource("Rcode/ReadAffyData.R");
		System.out.println("Executing script: " + scriptURL);
		exeR.executeScript(scriptURL);
		System.out.println("");
		
		String probe3 = exeR.getStringValue("probes[3]");
		System.out.println("Third sample: " + probe3);
		
		String ed4 = exeR.getStringValue("ed[4]");
		System.out.println("Fourth value: " + ed4);
		
		String sample2 = exeR.getStringValue("samp[2]");
		System.out.println("Second sample: " + sample2);
	}
}
