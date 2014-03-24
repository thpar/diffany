package be.svlandeg.diffany.r;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import org.rosuda.JRI.REXP;

/**
 * This class exploits an existing R bridge to process R commands from within Diffany.
 * The R commands can be issued in-line, or in a separate .R file, according to the specific definitions of the methods.
 * 
 * @author Sofie Van Landeghem
 */
public class ExecuteR
{

	private RBridge bridge;

	/**
	 * Create a new instance to execute R commando's, through a certain R bridge.
	 */
	public ExecuteR(RBridge bridge)
	{
		this.bridge = bridge;
	}

	/**
	 * Change the execution directory for the R commando's.
	 * @param dir_path
	 * @return the old path setting
	 */
	public String changeExecutionDir(String dir_path)
	{
		bridge.evaluate("path <- getwd()");
		String old_dir_path = getStringValue("path");
		
		// R wants to get forward slashed in the path 
		String new_dir_path = dir_path.replace("\\", "/");
		
		bridge.evaluate("setwd('" + new_dir_path + "')"); 
		
		return old_dir_path;
	}
	
	/**
	 * TODO documentation
	 */
	public String getStringValue(String variable)
	{
		REXP value = bridge.evaluate(variable);
		return value.asString();
	}
		
	/**
	 * TODO documentation
	 */
	public double getDoubleValue(String variable)
	{
		REXP value = bridge.evaluate(variable);
		return value.asDouble();
	}
	
	/**
	 * TODO documentation
	 * TODO: will this code work when packaged inside a jar or will we need to create a tmp file?
	 */
	public void executeScript(URL scriptURL)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(new File(scriptURL.toURI())));
			String line = reader.readLine();
			while (line != null)
			{
				bridge.evaluate(line);
				line = reader.readLine();
			}
			reader.close();
		}
		catch (IOException e)
		{
			System.out.println("Couldn't read R code : " + e.getMessage());
			return;
		}
		catch (URISyntaxException e)
		{
			System.out.println("Couldn't read R code : " + e.getMessage());
			return;
		}
	}

}
