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
		System.out.println(" Fetching ID mapping data: " + mappingURL);
		Map<String, Set<String>> idmapping = new MapID().getAllArrayMappings(new File(mappingURL.toURI())); 
		System.out.println("");
		
		String[] samples = exeR.getStringArray("samples");
		System.out.println("  Samples: " + samples.length);
		System.out.println("  e.g. second sample: " + samples[1]);
		System.out.println("");
		
		String[] probesets = exeR.getStringArray("probesets");
		System.out.println("  Probe sets: " + probesets.length);
		System.out.println("  e.g. third probe: set " + probesets[2] + " -> " + idmapping.get(probesets[2]));
		System.out.println("");
		
		double[][] expressionValues = exeR.getDoubleMatrix("expressionMatrix");
		System.out.println("  Expression values dimension: " + expressionValues.length + " - " + expressionValues[0].length);
		System.out.println("  e.g. third value in first two samples: " + expressionValues[2][0] + " & " + expressionValues[2][1]);
		System.out.println("");
	}
}
