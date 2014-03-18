package be.svlandeg.diffany.usecase.osmotic;

import java.io.File;
import java.io.IOException;

/**
 * This class reads and processes the raw input data.
 * 
 * @author Sofie Van Landeghem
 */
public class InputProcessing
{

	/**
	 * TODO documentation v2.0
	 */
	public void processOsmoticData(File osmoticStressDir) throws IOException
	{
		System.out.println("reading " + osmoticStressDir.getAbsolutePath());
		for (File f : osmoticStressDir.listFiles())
		{
			String fileName = f.getName();
			if (fileName.endsWith(".CEL"))
			{
				System.out.println(" processing " + fileName);
				
				// TODO: we'll use R and JRI to read and process the raw data !
			}
		}
		
	}

}
