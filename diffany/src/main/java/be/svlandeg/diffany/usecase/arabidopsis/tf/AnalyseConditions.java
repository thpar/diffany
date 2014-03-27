package be.svlandeg.diffany.usecase.arabidopsis.tf;

import java.util.Set;

import be.svlandeg.diffany.core.expression.ExpressionData;
import be.svlandeg.diffany.core.networks.GenericNetwork;

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
	public void integrateTFandExpr(GenericNetwork tfNetwork, Set<ExpressionData> datasets)
	{
		System.out.println(" Analysing " + tfNetwork.getStringRepresentation());
		System.out.println("  " + tfNetwork.getNodes().size() + " genes");
		System.out.println("  " + tfNetwork.getEdges().size() + " interactions");
		
		System.out.println("  ");
		
		System.out.println(" Analysing all expression datasets:");
		for (ExpressionData data : datasets)
		{
			System.out.println("  " + data);
		}
		
	}

}
