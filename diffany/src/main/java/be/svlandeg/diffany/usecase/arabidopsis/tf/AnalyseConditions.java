package be.svlandeg.diffany.usecase.arabidopsis.tf;

import java.util.Set;

import be.svlandeg.diffany.core.expression.ExpressionData;

/**
 * This class analyses the TF-target data and combines it with condition-specific expression data to determine
 * condition-dependent TF-target networks
 * 
 * @author Sofie Van Landeghem
 */
public class AnalyseConditions
{
	
	/**
	 * TODO v2.1 documentation
	 * @param datasets
	 */
	public void integrateTFandExpr(Set<ExpressionData> datasets)
	{
		System.out.println(" Analysing all expression datasets:");
		for (ExpressionData data : datasets)
		{
			System.out.println("  " + data);
		}
		
	}

}
