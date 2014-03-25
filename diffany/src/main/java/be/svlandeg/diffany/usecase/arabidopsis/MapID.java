package be.svlandeg.diffany.usecase.arabidopsis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;

/**
 * This class provides functionality to map Arabidopsis identifiers, for instance between array IDs and locus IDs
 * 
 * @author Sofie Van Landeghem
 */
public class MapID
{
	
	/**
	 * Obtain the mapping data between array IDs (keys) and locus IDs (values, as a set).
	 * 
	 * @param inputfile the .tab file containing array elements and a set of locus IDs, all separated by tabs.
	 * @return the mapping of array element IDs to their corresponding (set of) locus IDs
	 * @throws IOException when the input file can not be read properly
	 */
	public Map<String, Set<String>> getAllArrayMappings(File inputfile) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(inputfile));
		
		Map<String, Set<String>> map = new HashMap<String, Set<String>>();
		
		String line = reader.readLine();
		while (line != null)
		{
			StringTokenizer stok = new StringTokenizer(line, "\t");
			String array_element = stok.nextToken();
			Set<String> locusIDs = new HashSet<String>();
			
			if (map.containsKey(array_element))
			{
				System.out.println("error: found " + array_element + " twice");
			}
			
			while (stok.hasMoreTokens())
			{
				String locus = stok.nextToken();
				locusIDs.add(locus);
			}
			
			map.put(array_element, locusIDs);
			line = reader.readLine();
		}
		reader.close();
		return map;
	}

	/**
	 * Process the raw mapping data obtained from TAIR. 
	 * This should only be done once, from then on, the output file can be used for faster ID mapping.
	 * 
	 * @param inputfile the original .tab file containing array elements and a tab-delimited description, downloaded from TAIR
	 * @param outputfile the .tab file containing array elements and a set of locus IDs, all separated by tabs
	 * @throws IOException when the input file can not be read properly, or the output file can not be written
	 */
	protected void processRawATH1mapping(File inputfile, File outputfile) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputfile, false));
		
		BufferedReader reader = new BufferedReader(new FileReader(inputfile));
		String line = reader.readLine();
		while (line != null)
		{
			//System.out.println(line);
			StringTokenizer stok = new StringTokenizer(line, "\t");
			String array_element = stok.nextToken();
			stok.nextToken();					// type - information not used
			String organism = stok.nextToken();
			stok.nextToken(); 					// control - information not used
			String locus = stok.nextToken();
			if (organism.equals("Arabidopsis thaliana"))
			{
				StringTokenizer locusStok = new StringTokenizer(locus, ";");
				writer.append(array_element);
				writer.append("\t");
				while (locusStok.hasMoreTokens())
				{
					writer.append(locusStok.nextToken());
					writer.append("\t");
				}
				writer.newLine();
				writer.flush();
			}
			line = reader.readLine();
		}
		writer.flush();
		writer.close();
		
		reader.close();
	}

	/**
	 * Process the raw mapping data obtained from 
	 * ftp://ftp.arabidopsis.org/home/tair/Microarrays/Affymetrix/affy_ATH1_array_elements-2010-12-20.txt
	 * 
	 * 
	 * Currently, the output (written to the D disk), has been stored as a resources file in Maven.
	 * There is thus no need to execute this method again.
	 */
	public static void main(String[] args) throws IOException
	{
		String inputRoot = "D:" + File.separator + "diffany-osmotic";
		
		File inputfile = new File(inputRoot, "affy_ATH1_array_elements.tab");
		File outputfile = new File(inputRoot, "affy_ATH1_ID_mapping.tab");
		
		System.out.println("Transferring " + inputfile + " to " + outputfile);
		new MapID().processRawATH1mapping(inputfile, outputfile);
		System.out.println("Done !");
	}

}
