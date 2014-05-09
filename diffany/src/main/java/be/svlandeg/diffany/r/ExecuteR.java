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
	 * @param bridge the bridge that allows talking to R
	 */
	public ExecuteR(RBridge bridge)
	{
		this.bridge = bridge;
	}
	
	/**
	 * Change the execution directory for the R commando's.
	 * 
	 * @param dir_path the new working directory
	 * @return the old working directory
	 */
	public String changeExecutionDir(String dir_path)
	{
		bridge.evaluate("path <- getwd()");
		String old_dir_path = getStringValue("path");
		String new_dir_path = bridge.convertSlashes(dir_path);

		bridge.evaluate("setwd('" + new_dir_path + "')");

		return old_dir_path;
	}
	
	/**
	 * Retrieve whether a variable by the given name exists or not.
	 * 
	 * @param variable the symbol of the variable (as defined previously)
	 * @return whether or not the variable exists in the R environment.
	 */
	public boolean doesVariableExist(String variable)
	{
		REXP value = bridge.evaluate(variable);
		return value != null;
	}

	/**
	 * Retrieve the value of a certain String variable previously defined/calculated in R.
	 * 
	 * @param variable the symbol of the variable (as defined previously)
	 * @return the String value in the R environment
	 */
	public String getStringValue(String variable)
	{
		REXP value = bridge.evaluate(variable);
		if (doesVariableExist(variable))
		{
			return value.asString();
		}
		return null;
	}

	/**
	 * Retrieve the value of a certain String array previously defined/calculated in R.
	 * To check whether a variable exists, the method doesVariableExist is used.
	 * 
	 * @param variable the symbol of the variable (as defined previously)
	 * @return the String array in the R environment, or null if it doesn't exist in the environment
	 */
	public String[] getStringArray(String variable)
	{
		REXP value = bridge.evaluate(variable);
		if (doesVariableExist(variable))
		{
			return value.asStringArray();
		}
		return null;
	}

	/**
	 * Retrieve the value of a certain Double variable previously defined/calculated in R.
	 * To check whether a variable exists, the method doesVariableExist is used.
	 * 
	 * @param variable the symbol of the variable (as defined previously)
	 * @return the Double value in the R environment, or null if it doesn't exist in the environment
	 */
	public Double getDoubleValue(String variable)
	{
		REXP value = bridge.evaluate(variable);
		if (doesVariableExist(variable))
		{
			return value.asDouble();
		}
		return null;
	}

	/**
	 * Retrieve the value of a certain double array previously defined/calculated in R.
	 * To check whether a variable exists, the method doesVariableExist is used.
	 * 
	 * @param variable the symbol of the variable (as defined previously)
	 * @return the double array in the R environment, or null if it doesn't exist in the environment
	 */
	public double[] getDoubleArray(String variable)
	{
		REXP value = bridge.evaluate(variable);
		if (doesVariableExist(variable))
		{
			return value.asDoubleArray();
		}
		return null;
	}

	
	/**
	 * Retrieve the value of a certain double matrix previously defined/calculated in R.
	 * To check whether a variable exists, the method doesVariableExist is used.
	 * 
	 * @param variable the symbol of the variable (as defined previously)
	 * @return the double matrix in the R environment, or null if it doesn't exist in the environment
	 */
	public double[][] getDoubleMatrix(String variable)
	{
		REXP value = bridge.evaluate(variable);
		if (doesVariableExist(variable))
		{
			return value.asDoubleMatrix();
		}
		return null;
	}

	/**
	 * Execute a certain script in R by evaluating each line consequently.
	 * 
	 * @param scriptURL the URL (location) of the script that needs to be executed.
	 * @throws URISyntaxException
	 * @throws IOException
	 */
	public void executeScript(URL scriptURL) throws URISyntaxException, IOException
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

}
