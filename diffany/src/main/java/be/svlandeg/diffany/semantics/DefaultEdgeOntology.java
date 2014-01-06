package be.svlandeg.diffany.semantics;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;

import be.svlandeg.diffany.concepts.VisualEdgeStyle.ArrowHead;

/**
 * This class provides initial suggestions to the user concerning edge type-to-category mappings.
 * 
 * @author Sofie Van Landeghem
 */
public class DefaultEdgeOntology extends TreeEdgeOntology
{
	

	/**
	 * Create a default edge ontology, with default edge categories and type-category mappings.
	 */
	public DefaultEdgeOntology()
	{
		super ("increase", "increases", "decrease", "decreases", defineAllCategories());
		
		insertDefaultParents();
		insertDefaultMappings();
		insertDefaultTranslations();
		addDefaultPaintedParents();
	}
	
	/**
	 * Define all categories to be used in this EdgeOntology
	 * 
	 * @return the defined categories
	 */
	protected static Map<String, Boolean> defineAllCategories()
	{
		Map<String, Boolean> cats = new HashMap<String, Boolean>();
		
		cats.put("genetic_interaction", true);
		cats.put("positive_genetic_interaction", true);
		cats.put("negative_genetic_interaction", true);
		cats.put("synthetic_lethality", true);
		
		cats.put("regulation", false);
		cats.put("positive_regulation", false);
		cats.put("catalysis", false);
		cats.put("negative_regulation", false);
		cats.put("inhibition", false);

		cats.put("ppi", true);
		cats.put("colocalization", true);
		cats.put("coexpression", true);
		
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
	 * Provide some default child-parent mapping between source categories.
	 */
	protected void insertDefaultParents()
	{
		putSourceParent("positive_regulation", "regulation");
		putSourceParent("negative_regulation", "regulation");
		putSourceParent("catalysis", "positive_regulation");
		putSourceParent("inhibition", "negative_regulation");
		
		putSourceParent("positive_genetic_interaction", "genetic_interaction");
		putSourceParent("negative_genetic_interaction", "genetic_interaction");
		putSourceParent("synthetic_lethality", "negative_genetic_interaction");
		
		putSourceParent("phosphorylation", "ptm");
		putSourceParent("dephosphorylation", "ptm");
		putSourceParent("ubiquitination", "ptm");
		putSourceParent("deubiquitination", "ptm");
		putSourceParent("methylation", "ptm");
		putSourceParent("demethylation", "ptm");
		putSourceParent("hydroxylation", "ptm");
		putSourceParent("dehydroxylation", "ptm");
		putSourceParent("acetylation", "ptm");
		putSourceParent("deacetylation", "ptm");
		putSourceParent("glycosylation", "ptm");
		putSourceParent("deglycosylation", "ptm");
	}
	
	/**
	 * Provide some default translation rules.
	 * TODO: more elegant solution for this?!
	 */
	protected void insertDefaultTranslations()
	{
		/// REGULATION
		
		putDifferentialTranslation("positive_regulation_to_negative_regulation", "decreases_regulation");
		putDifferentialTranslation("positive_regulation_to_inhibition", "decreases_regulation");
		putDifferentialTranslation("catalysis_to_negative_regulation", "decreases_regulation");
		putDifferentialTranslation("catalysis_to_inhibition", "decreases_regulation");
		
		putDifferentialTranslation("negative_regulation_to_positive_regulation", "increases_regulation");
		putDifferentialTranslation("negative_regulation_to_catalysis", "increases_regulation");
		putDifferentialTranslation("inhibition_to_positive_regulation", "increases_regulation");
		putDifferentialTranslation("inhibition_to_catalysis", "increases_regulation");
		
		putDifferentialTranslation("decreases_positive_regulation", "decreases_regulation");
		putDifferentialTranslation("decreases_catalysis", "decreases_regulation");
		putDifferentialTranslation("increases_positive_regulation", "increases_regulation");
		putDifferentialTranslation("increases_catalysis", "increases_regulation");
		
		putDifferentialTranslation("decreases_negative_regulation", "increases_regulation");
		putDifferentialTranslation("decreases_inhibition", "increases_regulation");
		putDifferentialTranslation("increases_negative_regulation", "decreases_regulation");
		putDifferentialTranslation("increases_inhibition", "decreases_regulation");

		/// GI
		
		putDifferentialTranslation("positive_genetic_interaction_to_negative_genetic_interaction", "decrease_genetic_interaction");
		putDifferentialTranslation("positive_genetic_interaction_to_synthetic_lethality", "decrease_genetic_interaction");
		putDifferentialTranslation("negative_genetic_interaction_to_positive_genetic_interaction", "increase_genetic_interaction");
		putDifferentialTranslation("synthetic_lethality_to_positive_genetic_interaction", "increase_genetic_interaction");
		
		putDifferentialTranslation("decrease_positive_genetic_interaction", "decrease_genetic_interaction");
		putDifferentialTranslation("increase_positive_genetic_interaction", "increase_genetic_interaction");
		
		putDifferentialTranslation("decrease_negative_genetic_interaction", "increase_genetic_interaction");
		putDifferentialTranslation("increase_negative_genetic_interaction", "decrease_genetic_interaction");
		putDifferentialTranslation("decrease_synthetic_lethality", "increase_genetic_interaction");
		putDifferentialTranslation("increase_synthetic_lethality", "decrease_genetic_interaction");
	}
	
	/**
	 * Provide a default mapping edge type to category mapping.
	 */
	protected void addDefaultPaintedParents()
	{
		addColor("ptm", Color.BLUE);
		addColor("positive_regulation", Color.GREEN);
		addColor("negative_regulation", Color.RED);
		
		addArrowHead("ptm", ArrowHead.ARROW);
	}

	/**
	 * Provide a default mapping edge type to category mapping.
	 */
	protected void insertDefaultMappings()
	{
		boolean overwrite = false;

		// PTM
		
		addSourceCategoryMapping("methylation", "methylation", overwrite);
		addSourceCategoryMapping("methylates", "methylation", overwrite);
		addSourceCategoryMapping("methylate", "methylation", overwrite);
		
		addSourceCategoryMapping("demethylation", "demethylation", overwrite);
		addSourceCategoryMapping("demethylates", "demethylation", overwrite);
		addSourceCategoryMapping("demethylate", "demethylation", overwrite);
		
		addSourceCategoryMapping("phosphorylation", "phosphorylation", overwrite);
		addSourceCategoryMapping("phosphorylates", "phosphorylation", overwrite);
		addSourceCategoryMapping("phosphorylate", "phosphorylation", overwrite);
		
		addSourceCategoryMapping("dephosphorylation", "dephosphorylation", overwrite);
		addSourceCategoryMapping("dephosphorylates", "dephosphorylation", overwrite);
		addSourceCategoryMapping("dephosphorylate", "dephosphorylation", overwrite);
		
		addSourceCategoryMapping("ubiquitination", "ubiquitination", overwrite);
		addSourceCategoryMapping("ubiquitinates", "ubiquitination", overwrite);
		addSourceCategoryMapping("ubiquitinate", "ubiquitination", overwrite);
		
		addSourceCategoryMapping("deubiquitination", "deubiquitination", overwrite);
		addSourceCategoryMapping("deubiquitinates", "deubiquitination", overwrite);
		addSourceCategoryMapping("deubiquitinate", "deubiquitination", overwrite);
		
		addSourceCategoryMapping("deglycosylation", "deglycosylation", overwrite);
		addSourceCategoryMapping("deglycosylates", "deglycosylation", overwrite);
		addSourceCategoryMapping("deglycosylate", "deglycosylation", overwrite);
		
		addSourceCategoryMapping("glycosylation", "glycosylation", overwrite);
		addSourceCategoryMapping("glycosylates", "glycosylation", overwrite);
		addSourceCategoryMapping("glycosylate", "glycosylation", overwrite);
		
		addSourceCategoryMapping("deacetylation", "deacetylation", overwrite);
		addSourceCategoryMapping("deacetylates", "deacetylation", overwrite);
		addSourceCategoryMapping("deacetylate", "deacetylation", overwrite);
		
		addSourceCategoryMapping("acetylation", "acetylation", overwrite);
		addSourceCategoryMapping("acetylates", "acetylation", overwrite);
		addSourceCategoryMapping("acetylate", "acetylation", overwrite);
		
		addSourceCategoryMapping("hydroxylation", "hydroxylation", overwrite);
		addSourceCategoryMapping("hydroxylates", "hydroxylation", overwrite);
		addSourceCategoryMapping("hydroxylate", "hydroxylation", overwrite);
		
		addSourceCategoryMapping("dehydroxylation", "dehydroxylation", overwrite);
		addSourceCategoryMapping("dehydroxylates", "dehydroxylation", overwrite);
		addSourceCategoryMapping("dehydroxylate", "dehydroxylation", overwrite);
		
		addSourceCategoryMapping("ptm", "ptm", overwrite);
		addSourceCategoryMapping("post-translational modification", "ptm", overwrite);
		addSourceCategoryMapping("post-translational_modification", "ptm", overwrite);
		addSourceCategoryMapping("post_translational_modification", "ptm", overwrite);
		addSourceCategoryMapping("post_translational modification", "ptm", overwrite);
		addSourceCategoryMapping("post-translational-modification", "ptm", overwrite);
		addSourceCategoryMapping("post translational modification", "ptm", overwrite);
		addSourceCategoryMapping("posttranslational modification", "ptm", overwrite);
		addSourceCategoryMapping("posttranslationalmodification", "ptm", overwrite);

		// binding and other symmetricals
		
		addSourceCategoryMapping("complex", "ppi", overwrite);
		addSourceCategoryMapping("complex formation", "ppi", overwrite);
		addSourceCategoryMapping("complex-formation", "ppi", overwrite);
		addSourceCategoryMapping("complex_formation", "ppi", overwrite);
		addSourceCategoryMapping("complexformation", "ppi", overwrite);
		addSourceCategoryMapping("ppi", "ppi", overwrite);
		addSourceCategoryMapping("protein-protein interaction", "ppi", overwrite);
		addSourceCategoryMapping("proteinprotein interaction", "ppi", overwrite);
		addSourceCategoryMapping("proteinproteininteraction", "ppi", overwrite);
		addSourceCategoryMapping("binds", "ppi", overwrite);
		addSourceCategoryMapping("bind", "ppi", overwrite);
		addSourceCategoryMapping("binding", "ppi", overwrite);
		
		addSourceCategoryMapping("transcription", "transcription", overwrite);
		addSourceCategoryMapping("transcribes", "transcription", overwrite);
		addSourceCategoryMapping("transcribe", "transcription", overwrite);
		addSourceCategoryMapping("tf", "transcription", overwrite);
		addSourceCategoryMapping("transcription factor", "transcription", overwrite);
		addSourceCategoryMapping("transcription-factor", "transcription", overwrite);
		addSourceCategoryMapping("transcription_factor", "transcription", overwrite);
		addSourceCategoryMapping("transcriptionfactor", "transcription", overwrite);
		addSourceCategoryMapping("transcription factor binding", "transcription", overwrite);
		addSourceCategoryMapping("transcription-factor binding", "transcription", overwrite);
		addSourceCategoryMapping("transcription_factor binding", "transcription", overwrite);
		addSourceCategoryMapping("tf binding", "transcription", overwrite);
		addSourceCategoryMapping("tf-binding", "transcription", overwrite);
		addSourceCategoryMapping("tf_binding", "transcription", overwrite);
		
		addSourceCategoryMapping("coexpression", "coexpression", overwrite);
		addSourceCategoryMapping("coexpressed", "coexpression", overwrite);
		addSourceCategoryMapping("coexpresses", "coexpression", overwrite);
		addSourceCategoryMapping("coexpress", "coexpression", overwrite);
		addSourceCategoryMapping("coexpresses with", "coexpression", overwrite);
		
		addSourceCategoryMapping("colocalization", "colocalization", overwrite);
		addSourceCategoryMapping("colocalized", "colocalization", overwrite);
		addSourceCategoryMapping("colocalizes", "colocalization", overwrite);
		addSourceCategoryMapping("colocalize", "colocalization", overwrite);
		addSourceCategoryMapping("colocaliz with", "colocalization", overwrite);
		addSourceCategoryMapping("colocalisation", "colocalization", overwrite);
		addSourceCategoryMapping("colocalised", "colocalization", overwrite);
		addSourceCategoryMapping("colocalises", "colocalization", overwrite);
		addSourceCategoryMapping("colocalise", "colocalization", overwrite);
		addSourceCategoryMapping("colocalises with", "colocalization", overwrite);

		// regulation category and common synonyms
		addSourceCategoryMapping("regulation", "regulation", overwrite);
		addSourceCategoryMapping("regulates", "regulation", overwrite);
		addSourceCategoryMapping("regulate", "regulation", overwrite);
		addSourceCategoryMapping("influence", "regulation", overwrite);
		addSourceCategoryMapping("influences", "regulation", overwrite);
		addSourceCategoryMapping("effect", "regulation", overwrite);
		addSourceCategoryMapping("effects", "regulation", overwrite);

		// positive regulation category and common synonyms
		addSourceCategoryMapping("positive regulation", "positive_regulation", overwrite);
		addSourceCategoryMapping("positive_regulation", "positive_regulation", overwrite);
		addSourceCategoryMapping("positive-regulation", "positive_regulation", overwrite);
		addSourceCategoryMapping("positiveregulation", "positive_regulation", overwrite);
		addSourceCategoryMapping("positive reg", "positive_regulation", overwrite);
		addSourceCategoryMapping("positive-reg", "positive_regulation", overwrite);
		addSourceCategoryMapping("positive_reg", "positive_regulation", overwrite);
		addSourceCategoryMapping("positivereg", "positive_regulation", overwrite);
		addSourceCategoryMapping("pos regulation", "positive_regulation", overwrite);
		addSourceCategoryMapping("pos_regulation", "positive_regulation", overwrite);
		addSourceCategoryMapping("pos-regulation", "positive_regulation", overwrite);
		addSourceCategoryMapping("posregulation", "positive_regulation", overwrite);
		addSourceCategoryMapping("pos reg", "positive_regulation", overwrite);
		addSourceCategoryMapping("pos-reg", "positive_regulation", overwrite);
		addSourceCategoryMapping("pos_reg", "positive_regulation", overwrite);
		addSourceCategoryMapping("posreg", "positive_regulation", overwrite);
		addSourceCategoryMapping("pos", "positive_regulation", overwrite);
		addSourceCategoryMapping("positive", "positive_regulation", overwrite);
		addSourceCategoryMapping("positively regulates", "positive_regulation", overwrite);
		addSourceCategoryMapping("positively_regulates", "positive_regulation", overwrite);
		addSourceCategoryMapping("positively-regulates", "positive_regulation", overwrite);
		addSourceCategoryMapping("positivelyregulates", "positive_regulation", overwrite);
		addSourceCategoryMapping("upregulation", "positive_regulation", overwrite);
		addSourceCategoryMapping("up-regulation", "positive_regulation", overwrite);
		addSourceCategoryMapping("up regulation", "positive_regulation", overwrite);
		
		addSourceCategoryMapping("catalysis", "catalysis", overwrite);
		addSourceCategoryMapping("catalyses", "catalysis", overwrite);
		addSourceCategoryMapping("catalysation", "catalysis", overwrite);
		addSourceCategoryMapping("catalyzation", "catalysis", overwrite);

		// negative regulation category and common synonyms
		addSourceCategoryMapping("negative regulation", "negative_regulation", overwrite);
		addSourceCategoryMapping("negative_regulation", "negative_regulation", overwrite);
		addSourceCategoryMapping("negative-regulation", "negative_regulation", overwrite);
		addSourceCategoryMapping("negativeregulation", "negative_regulation", overwrite);
		addSourceCategoryMapping("negative reg", "negative_regulation", overwrite);
		addSourceCategoryMapping("negative_reg", "negative_regulation", overwrite);
		addSourceCategoryMapping("negative-reg", "negative_regulation", overwrite);
		addSourceCategoryMapping("negativereg", "negative_regulation", overwrite);
		addSourceCategoryMapping("neg regulation", "negative_regulation", overwrite);
		addSourceCategoryMapping("neg_regulation", "negative_regulation", overwrite);
		addSourceCategoryMapping("neg-regulation", "negative_regulation", overwrite);
		addSourceCategoryMapping("negregulation", "negative_regulation", overwrite);
		addSourceCategoryMapping("neg reg", "negative_regulation", overwrite);
		addSourceCategoryMapping("neg_reg", "negative_regulation", overwrite);
		addSourceCategoryMapping("neg-reg", "negative_regulation", overwrite);
		addSourceCategoryMapping("negreg", "negative_regulation", overwrite);
		addSourceCategoryMapping("neg", "negative_regulation", overwrite);
		addSourceCategoryMapping("negative", "negative_regulation", overwrite);
		addSourceCategoryMapping("negatively regulates", "negative_regulation", overwrite);
		addSourceCategoryMapping("negatively_regulates", "negative_regulation", overwrite);
		addSourceCategoryMapping("negatively-regulates", "negative_regulation", overwrite);
		addSourceCategoryMapping("negativelyregulates", "negative_regulation", overwrite);
		addSourceCategoryMapping("downregulation", "negative_regulation", overwrite);
		addSourceCategoryMapping("down-regulation", "negative_regulation", overwrite);
		addSourceCategoryMapping("down regulation", "negative_regulation", overwrite);
		
		addSourceCategoryMapping("inhibits", "inhibition", overwrite);
		addSourceCategoryMapping("inhibit", "inhibition", overwrite);
		addSourceCategoryMapping("inhibition", "inhibition", overwrite);
		
		// genetic interactions
		addSourceCategoryMapping("positive_genetic_interaction", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("positive genetic interaction", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("positive-genetic-interaction", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("positive gi", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("alleviating gi", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("alleviating genetic_interaction", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("alleviating_genetic_interaction", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("alleviating-genetic-interaction", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("positive_epistatic_interaction", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("positive epistatic interaction", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("positive-epistatic-interaction", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("alleviating_epistatic_interaction", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("alleviating epistatic interaction", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("alleviating-epistatic-interaction", "positive_genetic_interaction", overwrite);
		
		addSourceCategoryMapping("negative_genetic_interaction", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("negative genetic interaction", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("negative-genetic-interaction", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("negative gi", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("aggrevating_genetic_interaction", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("aggrevating genetic interaction", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("aggrevating-genetic-interaction", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("aggrevating gi", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("negative_epistatic_interaction", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("negative epistatic interaction", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("negative-epistatic-interaction", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("aggrevating_epistatic_interaction", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("aggrevating epistatic interaction", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("aggrevating-epistatic-interaction", "negative_genetic_interaction", overwrite);
		
		addSourceCategoryMapping("synthetic_lethality", "synthetic_lethality", overwrite);
		addSourceCategoryMapping("synthetic lethality", "synthetic_lethality", overwrite);
		addSourceCategoryMapping("synthetic-lethality", "synthetic_lethality", overwrite);
		addSourceCategoryMapping("synthetically_lethal", "synthetic_lethality", overwrite);
		addSourceCategoryMapping("synthetically lethal", "synthetic_lethality", overwrite);
		addSourceCategoryMapping("sl", "synthetic_lethality", overwrite);
		
		addSourceCategoryMapping("genetic_interaction", "genetic_interaction", overwrite);
		addSourceCategoryMapping("gi", "genetic_interaction", overwrite);
		addSourceCategoryMapping("epistatic", "genetic_interaction", overwrite);
		addSourceCategoryMapping("epistasis", "genetic_interaction", overwrite);
		addSourceCategoryMapping("epistatic interaction", "genetic_interaction", overwrite);
	}

}
