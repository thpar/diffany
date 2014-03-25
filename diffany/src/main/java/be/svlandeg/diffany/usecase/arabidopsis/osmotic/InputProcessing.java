package be.svlandeg.diffany.usecase.arabidopsis.osmotic;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.r.ExecuteR;
import be.svlandeg.diffany.usecase.arabidopsis.MapID;

/**
 * This class reads and processes the raw input data.
 * 
 * @author Sofie Van Landeghem
 */
public class InputProcessing
{

	/**
	 * Process the raw osmotic expression data with R.
	 * 
	 * Currently, the needed R script is loaded from the context, and is defined in the 'resources' folder
	 * of the Maven project.
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
		
		// TODO V2.0: currently this assumes libs "affy", "affyPLM" and "org.Dm.eg.db" are pre-installed!
		URL scriptURL = Thread.currentThread().getContextClassLoader().getResource("Rcode/ReadAffyData.R");
		System.out.println(" Executing script: " + scriptURL);
		System.out.println(" (this may take a minute ... please be patient and do not interrupt the execution) ");
		exeR.executeScript(scriptURL);
		System.out.println("");
		
		URL mappingURL = Thread.currentThread().getContextClassLoader().getResource("data/affy_ATH1_ID_mapping.tab");
		Map<String, Set<String>> idmapping = new MapID().getAllArrayMappings(new File(mappingURL.toURI())); 
		
		String[] probes = exeR.getStringArray("probes");
		System.out.println("  Third probe: " + probes[2] + " -> " + idmapping.get(probes[2]));
		
		double[] values = exeR.getDoubleArray("values");
		System.out.println("  Fourth value: " + values[3]);
		
		String[] samples = exeR.getStringArray("samples");
		System.out.println("  Second sample: " + samples[1]);
	}
}
