package be.svlandeg.diffany.examples;

import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

public class OsmoticUseCase extends GenericExample{

	private NodeMapper nm;

	public OsmoticUseCase() {
		nm = new DefaultNodeMapper();
	}
	
	public Project getTestProject(){
		String name = "Osmotic Use Case";
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo, nm);
		return p;
	}
	
	public int getTestConfiguration(Project p, int supportingCutoff){
		ReferenceNetwork r = getTestReference();
		Set<ConditionNetwork> c = getTestConditions();
		boolean cleanInput = true;
		int ID = p.addRunConfiguration(r,  c, supportingCutoff, cleanInput, null);
		return ID;
	}

	private ReferenceNetwork getTestReference() {
		// TODO Auto-generated method stub
		return null;
	}
	
	private Set<ConditionNetwork> getTestConditions() {
		Set<ConditionNetwork> cnetworks = new HashSet<ConditionNetwork>();
		cnetworks.add(getFirstCondition());
		cnetworks.add(getSecondCondition());
		cnetworks.add(getThirdCondition());
		cnetworks.add(getFourthCondition());
		return cnetworks;
	}

	private ConditionNetwork getFirstCondition() {
		String description = "1.5h";
		return null;
	}
	private ConditionNetwork getSecondCondition() {
		String description = "3h";
		return null;
	}
	private ConditionNetwork getThirdCondition() {
		String description = "12h";
		return null;
	}	
	private ConditionNetwork getFourthCondition() {
		String description = "24h";
		return null;
	}




	
}
