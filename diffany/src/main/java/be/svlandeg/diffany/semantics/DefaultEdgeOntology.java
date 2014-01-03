package be.svlandeg.diffany.semantics;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.concepts.EdgeDefinition;
import be.svlandeg.diffany.concepts.VisualEdgeStyle.ArrowHead;

/**
 * This class holds both an ActivityFlowEdgeOntology and a ProcessEdgeOntology. 
 * It can be used to give initial suggestions to the user concerning edge type-to-category mappings.
 * 
 * @author Sofie Van Landeghem
 */
public class DefaultEdgeOntology extends EdgeOntology
{

	private ActivityFlowEdgeOntology afOntology;
	private ProcessEdgeOntology prOntology;

	/**
	 * Create a default edge ontology, with generic up/down edge categories and
	 * type-category mappings.
	 */
	public DefaultEdgeOntology()
	{
		String pos_cat_symm = "increase";
		String pos_cat_dir = "increases";
		Map<String, Boolean> pos_cats = definePosCategories();

		String neg_cat_symm = "decrease";
		String neg_cat_dir = "decreases";
		Map<String, Boolean> neg_cats = defineNegCategories();
		
		Map<String, Boolean> neutral_cats = defineNeutralCategories();
		
		afOntology = new ActivityFlowEdgeOntology(pos_cat_symm, pos_cat_dir, neg_cat_symm, neg_cat_dir, pos_cats, neg_cats, neutral_cats);

		Map<String, Boolean> process_cats = defineProcessCategories();

		prOntology = new ProcessEdgeOntology(pos_cat_symm + "_", pos_cat_dir + "_", neg_cat_symm + "_", neg_cat_dir + "_", process_cats);

		insertDefaultParents();
		insertDefaultMappings();
		addDefaultPaintedParents();
	}

	
	
	/**
	 * Get the semantic category of a certain edge type. Matching is done independent of upper/lower casing.
	 * Both the internal ActivityFlowEdgeOntology as well as the ProcessEdgeOntology are consulted for the mapping.
	 * 
	 * @param edgeType the original type of the edge in a network
	 * @return the semantic category of that edge or null if it is not mapped in this ontology
	 */
	public String getSourceCategory(String edgeType)
	{
		String hit = afOntology.getSourceCategory(edgeType);
		if (hit == null)
		{
			hit = prOntology.getSourceCategory(edgeType);
		}
		return hit;
	}
	
	/**
	 * Return all categories present in this ontology, 
	 * consulting both the internal ActivityFlowEdgeOntology as well as the ProcessEdgeOntology.
	 * @return the set of categories mapped in this ontology.
	 */
	public Set<String> getAllSourceCategories()
	{
		Set<String> allCats = afOntology.getAllSourceCategories();
		allCats.addAll(prOntology.getAllSourceCategories());
		return allCats;
	}
	
	/**
	 * Return all categories present in this ontology, 
	 * consulting both the internal ActivityFlowEdgeOntology as well as the ProcessEdgeOntology.
	 * @return the set of categories mapped in this ontology.
	 */
	public Set<String> getAllDiffCategories()
	{
		Set<String> allCats = afOntology.getAllDiffCategories();
		allCats.addAll(prOntology.getAllDiffCategories());
		return allCats;
	}

	@Override
	public EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, Set<EdgeDefinition> conEdges, double cutoff) throws IllegalArgumentException
	{
		String refType = refEdge.getType();
		
		int conCount = conEdges.size();
		int afCount = 0;
		int prCount = 0;
		
		for (EdgeDefinition conEdge : conEdges)
		{
			String conType = conEdge.getType();
			if (afOntology.isDefinedSourceType(conType))
				afCount++;
			if (prOntology.isDefinedSourceType(conType))
				prCount++;
		}
		
		if (afOntology.isDefinedSourceType(refType) && afCount == conCount)
		{
			return afOntology.getDifferentialEdge(refEdge, conEdges, cutoff);
		}
		if (prOntology.isDefinedSourceType(refType) && prCount == conCount)
		{
			return prOntology.getDifferentialEdge(refEdge, conEdges, cutoff);
		}
		String errormsg = "The types '" + refType + "' and '" + conEdges.size() + "' conditional types could not be compared in this edge ontology!";
		throw new IllegalArgumentException(errormsg);
	}

	@Override
	public EdgeDefinition getOverlapEdge(Set<EdgeDefinition> edges, double cutoff, boolean minOperator) throws IllegalArgumentException
	{
		int conCount = edges.size();
		int afCount = 0;
		int prCount = 0;
		
		for (EdgeDefinition conEdge : edges)
		{
			String conType = conEdge.getType();
			if (afOntology.isDefinedSourceType(conType))
				afCount++;
			if (prOntology.isDefinedSourceType(conType))
				prCount++;
		}
		
		if (afCount == conCount)
		{
			return afOntology.getOverlapEdge(edges, cutoff, minOperator);
		}
		if (prCount == conCount)
		{
			return prOntology.getOverlapEdge(edges, cutoff, minOperator);
		}
		String errormsg = "The edge types could not be compared in this edge ontology!";
		throw new IllegalArgumentException(errormsg);
	}
	
	@Override
	public boolean isSymmetricalSourceType(String edgeType)
	{
		if (afOntology.isDefinedSourceType(edgeType))
			return afOntology.isSymmetricalSourceType(edgeType);
		return prOntology.isSymmetricalSourceType(edgeType);
	}

	@Override
	protected Color getDifferentialEdgeColor(String category)
	{
		if (afOntology.isDefinedDiffCategory(category))
			return afOntology.getDifferentialEdgeColor(category);
		return prOntology.getDifferentialEdgeColor(category);
	}

	@Override
	protected Color getSourceEdgeColor(String edgeType)
	{
		if (afOntology.isDefinedSourceType(edgeType))
			return afOntology.getSourceEdgeColor(edgeType);
		return prOntology.getSourceEdgeColor(edgeType);
	}
	
	@Override
	protected ArrowHead getDifferentialEdgeArrowHead(String category)
	{
		if (afOntology.isDefinedDiffCategory(category))
			return afOntology.getDifferentialEdgeArrowHead(category);
		return prOntology.getDifferentialEdgeArrowHead(category);
	}

	@Override
	protected ArrowHead getSourceEdgeArrowHead(String edgeType)
	{
		if (afOntology.isDefinedSourceType(edgeType))
			return afOntology.getSourceEdgeArrowHead(edgeType);
		return prOntology.getSourceEdgeArrowHead(edgeType);
	}
	
	/////////////////////////////// BELOW START THE DEFAULT DEFINITIONS ////////////////////
	
	/**
	 * Define the process categories to be used in the ProcessEdgeOntology
	 * 
	 * @return the defined process categories
	 */
	protected Map<String, Boolean> defineProcessCategories()
	{
		Map<String, Boolean> cats = new HashMap<String, Boolean>();
		
		cats.put("ppi", true);
		cats.put("colocalization", true);
		cats.put("coexpression", true);
		
		cats.put("transcription", false);
		
		cats.put("ptm", false);
		cats.put("phosphorylation", false);
		cats.put("dephosphorylation", false);
		cats.put("glycosylation", false);
		cats.put("deglycosylation", false);
		cats.put("acetylation", false);
		cats.put("deacetylation", false);
		cats.put("hydroxylation", false);
		cats.put("dehydroxylation", false);
		cats.put("ubiquitination", false);
		cats.put("deubiquitination", false);
		cats.put("methylation", false);
		cats.put("demethylation", false);
		
		return cats;
	}

	/**
	 * Define the positive categories to be used in the ActivityFlowEdgeOntology
	 * 
	 * @return the defined positive categories
	 */
	protected Map<String, Boolean> definePosCategories()
	{
		Map<String, Boolean> pos_cats = new HashMap<String, Boolean>();
		pos_cats.put("positive_genetic_interaction", true);
		
		pos_cats.put("positive_regulation", false);
		pos_cats.put("catalysis", false);
		return pos_cats;
	}

	/**
	 * Define the negative categories to be used in the ActivityFlowEdgeOntology
	 * 
	 * @return the defined negative categories
	 */
	protected Map<String, Boolean> defineNegCategories()
	{
		Map<String, Boolean> neg_cats = new HashMap<String, Boolean>();
		neg_cats.put("negative_genetic_interaction", true);
		neg_cats.put("synthetic_lethality", true);
		
		neg_cats.put("negative_regulation", false);
		neg_cats.put("inhibition", false);
		return neg_cats;
	}
	
	/**
	 * Define the neutral categories to be used in the ActivityFlowEdgeOntology
	 * 
	 * @return the defined negative categories
	 */
	protected Map<String, Boolean> defineNeutralCategories()
	{
		Map<String, Boolean> neutral_cats = new HashMap<String, Boolean>();
		neutral_cats.put("genetic_interaction", true);
		neutral_cats.put("regulation", false);
		return neutral_cats;
	}
	
	/**
	 * Provide a default child-parent mapping between source categories.
	 */
	protected void insertDefaultParents()
	{
		afOntology.putSourceParent("positive_regulation", "regulation");
		afOntology.putSourceParent("negative_regulation", "regulation");
		afOntology.putSourceParent("catalysis", "positive_regulation");
		afOntology.putSourceParent("inhibition", "negative_regulation");
		
		afOntology.putSourceParent("positive_genetic_interaction", "genetic_interaction");
		afOntology.putSourceParent("negative_genetic_interaction", "genetic_interaction");
		afOntology.putSourceParent("synthetic_lethality", "negative_genetic_interaction");
		
		prOntology.putSourceParent("phosphorylation", "ptm");
		prOntology.putSourceParent("dephosphorylation", "ptm");
		prOntology.putSourceParent("ubiquitination", "ptm");
		prOntology.putSourceParent("deubiquitination", "ptm");
		prOntology.putSourceParent("methylation", "ptm");
		prOntology.putSourceParent("demethylation", "ptm");
		prOntology.putSourceParent("hydroxylation", "ptm");
		prOntology.putSourceParent("dehydroxylation", "ptm");
		prOntology.putSourceParent("acetylation", "ptm");
		prOntology.putSourceParent("deacetylation", "ptm");
		prOntology.putSourceParent("glycosylation", "ptm");
		prOntology.putSourceParent("deglycosylation", "ptm");
	}
	
	/**
	 * Provide a default mapping edge type to category mapping.
	 */
	protected void addDefaultPaintedParents()
	{
		prOntology.addColor("ptm", Color.BLUE);
		prOntology.addArrowHead("ptm", ArrowHead.ARROW);
	}

	/**
	 * Provide a default mapping edge type to category mapping.
	 */
	protected void insertDefaultMappings()
	{
		boolean overwrite = false;

		// PR categories and common synonyms: PTM
		
		prOntology.addSourceCategoryMapping("methylation", "methylation", overwrite);
		prOntology.addSourceCategoryMapping("methylates", "methylation", overwrite);
		prOntology.addSourceCategoryMapping("methylate", "methylation", overwrite);
		
		prOntology.addSourceCategoryMapping("demethylation", "demethylation", overwrite);
		prOntology.addSourceCategoryMapping("demethylates", "demethylation", overwrite);
		prOntology.addSourceCategoryMapping("demethylate", "demethylation", overwrite);
		
		prOntology.addSourceCategoryMapping("phosphorylation", "phosphorylation", overwrite);
		prOntology.addSourceCategoryMapping("phosphorylates", "phosphorylation", overwrite);
		prOntology.addSourceCategoryMapping("phosphorylate", "phosphorylation", overwrite);
		
		prOntology.addSourceCategoryMapping("dephosphorylation", "dephosphorylation", overwrite);
		prOntology.addSourceCategoryMapping("dephosphorylates", "dephosphorylation", overwrite);
		prOntology.addSourceCategoryMapping("dephosphorylate", "dephosphorylation", overwrite);
		
		prOntology.addSourceCategoryMapping("ubiquitination", "ubiquitination", overwrite);
		prOntology.addSourceCategoryMapping("ubiquitinates", "ubiquitination", overwrite);
		prOntology.addSourceCategoryMapping("ubiquitinate", "ubiquitination", overwrite);
		
		prOntology.addSourceCategoryMapping("deubiquitination", "deubiquitination", overwrite);
		prOntology.addSourceCategoryMapping("deubiquitinates", "deubiquitination", overwrite);
		prOntology.addSourceCategoryMapping("deubiquitinate", "deubiquitination", overwrite);
		
		prOntology.addSourceCategoryMapping("deglycosylation", "deglycosylation", overwrite);
		prOntology.addSourceCategoryMapping("deglycosylates", "deglycosylation", overwrite);
		prOntology.addSourceCategoryMapping("deglycosylate", "deglycosylation", overwrite);
		
		prOntology.addSourceCategoryMapping("glycosylation", "glycosylation", overwrite);
		prOntology.addSourceCategoryMapping("glycosylates", "glycosylation", overwrite);
		prOntology.addSourceCategoryMapping("glycosylate", "glycosylation", overwrite);
		
		prOntology.addSourceCategoryMapping("deacetylation", "deacetylation", overwrite);
		prOntology.addSourceCategoryMapping("deacetylates", "deacetylation", overwrite);
		prOntology.addSourceCategoryMapping("deacetylate", "deacetylation", overwrite);
		
		prOntology.addSourceCategoryMapping("acetylation", "acetylation", overwrite);
		prOntology.addSourceCategoryMapping("acetylates", "acetylation", overwrite);
		prOntology.addSourceCategoryMapping("acetylate", "acetylation", overwrite);
		
		prOntology.addSourceCategoryMapping("hydroxylation", "hydroxylation", overwrite);
		prOntology.addSourceCategoryMapping("hydroxylates", "hydroxylation", overwrite);
		prOntology.addSourceCategoryMapping("hydroxylate", "hydroxylation", overwrite);
		
		prOntology.addSourceCategoryMapping("dehydroxylation", "dehydroxylation", overwrite);
		prOntology.addSourceCategoryMapping("dehydroxylates", "dehydroxylation", overwrite);
		prOntology.addSourceCategoryMapping("dehydroxylate", "dehydroxylation", overwrite);
		
		prOntology.addSourceCategoryMapping("ptm", "ptm", overwrite);
		prOntology.addSourceCategoryMapping("post-translational modification", "ptm", overwrite);
		prOntology.addSourceCategoryMapping("post-translational_modification", "ptm", overwrite);
		prOntology.addSourceCategoryMapping("post_translational_modification", "ptm", overwrite);
		prOntology.addSourceCategoryMapping("post_translational modification", "ptm", overwrite);
		prOntology.addSourceCategoryMapping("post-translational-modification", "ptm", overwrite);
		prOntology.addSourceCategoryMapping("post translational modification", "ptm", overwrite);
		prOntology.addSourceCategoryMapping("posttranslational modification", "ptm", overwrite);
		prOntology.addSourceCategoryMapping("posttranslationalmodification", "ptm", overwrite);

		// PR categories and common synonyms: binding and other symmetricals
		
		prOntology.addSourceCategoryMapping("complex", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("complex formation", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("complex-formation", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("complex_formation", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("complexformation", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("ppi", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("protein-protein interaction", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("proteinprotein interaction", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("proteinproteininteraction", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("binds", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("bind", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("binding", "ppi", overwrite);
		
		prOntology.addSourceCategoryMapping("transcription", "transcription", overwrite);
		prOntology.addSourceCategoryMapping("transcribes", "transcription", overwrite);
		prOntology.addSourceCategoryMapping("transcribe", "transcription", overwrite);
		prOntology.addSourceCategoryMapping("tf", "transcription", overwrite);
		prOntology.addSourceCategoryMapping("transcription factor", "transcription", overwrite);
		prOntology.addSourceCategoryMapping("transcription-factor", "transcription", overwrite);
		prOntology.addSourceCategoryMapping("transcription_factor", "transcription", overwrite);
		prOntology.addSourceCategoryMapping("transcriptionfactor", "transcription", overwrite);
		prOntology.addSourceCategoryMapping("transcription factor binding", "transcription", overwrite);
		prOntology.addSourceCategoryMapping("transcription-factor binding", "transcription", overwrite);
		prOntology.addSourceCategoryMapping("transcription_factor binding", "transcription", overwrite);
		prOntology.addSourceCategoryMapping("tf binding", "transcription", overwrite);
		prOntology.addSourceCategoryMapping("tf-binding", "transcription", overwrite);
		prOntology.addSourceCategoryMapping("tf_binding", "transcription", overwrite);
		
		prOntology.addSourceCategoryMapping("coexpression", "coexpression", overwrite);
		prOntology.addSourceCategoryMapping("coexpressed", "coexpression", overwrite);
		prOntology.addSourceCategoryMapping("coexpresses", "coexpression", overwrite);
		prOntology.addSourceCategoryMapping("coexpress", "coexpression", overwrite);
		prOntology.addSourceCategoryMapping("coexpresses with", "coexpression", overwrite);
		
		prOntology.addSourceCategoryMapping("colocalization", "colocalization", overwrite);
		prOntology.addSourceCategoryMapping("colocalized", "colocalization", overwrite);
		prOntology.addSourceCategoryMapping("colocalizes", "colocalization", overwrite);
		prOntology.addSourceCategoryMapping("colocalize", "colocalization", overwrite);
		prOntology.addSourceCategoryMapping("colocaliz with", "colocalization", overwrite);
		prOntology.addSourceCategoryMapping("colocalisation", "colocalization", overwrite);
		prOntology.addSourceCategoryMapping("colocalised", "colocalization", overwrite);
		prOntology.addSourceCategoryMapping("colocalises", "colocalization", overwrite);
		prOntology.addSourceCategoryMapping("colocalise", "colocalization", overwrite);
		prOntology.addSourceCategoryMapping("colocalises with", "colocalization", overwrite);

		// regulation category and common synonyms
		afOntology.addSourceCategoryMapping("regulation", "regulation", overwrite);
		afOntology.addSourceCategoryMapping("regulates", "regulation", overwrite);
		afOntology.addSourceCategoryMapping("regulate", "regulation", overwrite);
		afOntology.addSourceCategoryMapping("influence", "regulation", overwrite);
		afOntology.addSourceCategoryMapping("influences", "regulation", overwrite);
		afOntology.addSourceCategoryMapping("effect", "regulation", overwrite);
		afOntology.addSourceCategoryMapping("effects", "regulation", overwrite);

		// positive regulation category and common synonyms
		afOntology.addSourceCategoryMapping("positive regulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positive_regulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positive-regulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positiveregulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positive reg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positive-reg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positive_reg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positivereg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("pos regulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("pos_regulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("pos-regulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("posregulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("pos reg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("pos-reg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("pos_reg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("posreg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("pos", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positive", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positively regulates", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positively_regulates", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positively-regulates", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positivelyregulates", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("upregulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("up-regulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("up regulation", "positive_regulation", overwrite);
		
		afOntology.addSourceCategoryMapping("catalysis", "catalysis", overwrite);
		afOntology.addSourceCategoryMapping("catalyses", "catalysis", overwrite);
		afOntology.addSourceCategoryMapping("catalysation", "catalysis", overwrite);
		afOntology.addSourceCategoryMapping("catalyzation", "catalysis", overwrite);

		// negative regulation category and common synonyms
		afOntology.addSourceCategoryMapping("negative regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negative_regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negative-regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negativeregulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negative reg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negative_reg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negative-reg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negativereg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("neg regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("neg_regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("neg-regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negregulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("neg reg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("neg_reg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("neg-reg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negreg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("neg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negative", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negatively regulates", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negatively_regulates", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negatively-regulates", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negativelyregulates", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("downregulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("down-regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("down regulation", "negative_regulation", overwrite);
		
		afOntology.addSourceCategoryMapping("inhibits", "inhibition", overwrite);
		afOntology.addSourceCategoryMapping("inhibit", "inhibition", overwrite);
		afOntology.addSourceCategoryMapping("inhibition", "inhibition", overwrite);
		
		// genetic interactions
		afOntology.addSourceCategoryMapping("positive_genetic_interaction", "positive_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("positive genetic interaction", "positive_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("positive-genetic-interaction", "positive_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("positive gi", "positive_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("alleviating gi", "positive_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("alleviating genetic_interaction", "positive_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("alleviating_genetic_interaction", "positive_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("alleviating-genetic-interaction", "positive_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("positive_epistatic_interaction", "positive_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("positive epistatic interaction", "positive_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("positive-epistatic-interaction", "positive_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("alleviating_epistatic_interaction", "positive_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("alleviating epistatic interaction", "positive_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("alleviating-epistatic-interaction", "positive_genetic_interaction", overwrite);
		
		afOntology.addSourceCategoryMapping("negative_genetic_interaction", "negative_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("negative genetic interaction", "negative_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("negative-genetic-interaction", "negative_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("negative gi", "negative_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("aggrevating_genetic_interaction", "negative_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("aggrevating genetic interaction", "negative_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("aggrevating-genetic-interaction", "negative_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("aggrevating gi", "negative_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("negative_epistatic_interaction", "negative_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("negative epistatic interaction", "negative_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("negative-epistatic-interaction", "negative_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("aggrevating_epistatic_interaction", "negative_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("aggrevating epistatic interaction", "negative_genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("aggrevating-epistatic-interaction", "negative_genetic_interaction", overwrite);
		
		afOntology.addSourceCategoryMapping("synthetic_lethality", "synthetic_lethality", overwrite);
		afOntology.addSourceCategoryMapping("synthetic lethality", "synthetic_lethality", overwrite);
		afOntology.addSourceCategoryMapping("synthetic-lethality", "synthetic_lethality", overwrite);
		afOntology.addSourceCategoryMapping("synthetically_lethal", "synthetic_lethality", overwrite);
		afOntology.addSourceCategoryMapping("synthetically lethal", "synthetic_lethality", overwrite);
		afOntology.addSourceCategoryMapping("sl", "synthetic_lethality", overwrite);
		
		afOntology.addSourceCategoryMapping("genetic_interaction", "genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("gi", "genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("epistatic", "genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("epistasis", "genetic_interaction", overwrite);
		afOntology.addSourceCategoryMapping("epistatic interaction", "genetic_interaction", overwrite);
	}

}
