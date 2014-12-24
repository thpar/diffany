package be.svlandeg.diffany.core.semantics;

/*
 * #%L
 * Diffany
 * %%
 * Copyright (C) 2014 PSB/UGent - Sofie Van Landeghem and Thomas Van Parys
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU General Lesser Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/lgpl-3.0.html>.
 * #L%
 */

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;

import be.svlandeg.diffany.core.visualstyle.DefaultDiffEdgeDrawing;
import be.svlandeg.diffany.core.visualstyle.DefaultSourceEdgeDrawing;
import be.svlandeg.diffany.core.visualstyle.EdgeDrawing;
import be.svlandeg.diffany.core.visualstyle.EdgeStyle.ArrowHead;

/**
 * This class provides initial suggestions to the user concerning edge type-to-category mappings.
 * 
 * @author Sofie Van Landeghem
 */
public class DefaultEdgeOntology extends TreeEdgeOntology
{

	protected DefaultSourceEdgeDrawing sourceDraw;
	protected DefaultDiffEdgeDrawing diffDraw;

	/**
	 * Create a default edge ontology, with default edge categories and type-category mappings.
	 * For drawing the edges, a {@link DefaultSourceEdgeDrawing} and a {@link DefaultDiffEdgeDrawing} are created.
	 */
	public DefaultEdgeOntology()
	{
		super("increase_", "increases_", "decrease_", "decreases_", "unspecified_", defineAllCategories());
		sourceDraw = new DefaultSourceEdgeDrawing(this);
		diffDraw = new DefaultDiffEdgeDrawing(this);

		insertDefaultParents();
		insertDefaultMappings();
		insertDefaultContrasts();
		addDefaultPaintedParents();
	}

	@Override
	public EdgeDrawing getDifferentialEdgeDrawing()
	{
		return diffDraw;
	}

	@Override
	public EdgeDrawing getSourceEdgeDrawing()
	{
		return sourceDraw;
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
		
		cats.put("expression", false);
		cats.put("underexpression", false);
		cats.put("overexpression", false);
		
		cats.put("coexpression", true);

		cats.put("ppi", true);
		cats.put("colocalization", true);

		cats.put("protein-dna_binding", false);
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
	 * All default root categories should be children of either GENERIC_DIRECTED_CAT or GENERIC_SYMMETRICAL_CAT
	 */
	protected void insertDefaultParents()
	{
		putSourceCatParent("regulation", GENERIC_DIRECTED_CAT);
		putSourceCatParent("positive_regulation", "regulation");
		putSourceCatParent("negative_regulation", "regulation");
		putSourceCatParent("catalysis", "positive_regulation");
		putSourceCatParent("inhibition", "negative_regulation");
		
		putSourceCatParent("expression", GENERIC_DIRECTED_CAT);
		putSourceCatParent("overexpression", "expression");
		putSourceCatParent("underexpression", "expression");

		putSourceCatParent("genetic_interaction", GENERIC_SYMMETRICAL_CAT);
		putSourceCatParent("positive_genetic_interaction", "genetic_interaction");
		putSourceCatParent("negative_genetic_interaction", "genetic_interaction");
		putSourceCatParent("synthetic_lethality", "negative_genetic_interaction");
		
		putSourceCatParent("coexpression", GENERIC_SYMMETRICAL_CAT);

		putSourceCatParent("colocalization", GENERIC_SYMMETRICAL_CAT);
		putSourceCatParent("ppi", "colocalization");
		
		putSourceCatParent("protein-dna_binding", GENERIC_DIRECTED_CAT);
		putSourceCatParent("transcription", "protein-dna_binding");

		putSourceCatParent("ptm", GENERIC_DIRECTED_CAT);
		putSourceCatParent("phosphorylation", "ptm");
		putSourceCatParent("dephosphorylation", "ptm");
		putSourceCatParent("ubiquitination", "ptm");
		putSourceCatParent("deubiquitination", "ptm");
		putSourceCatParent("methylation", "ptm");
		putSourceCatParent("demethylation", "ptm");
		putSourceCatParent("hydroxylation", "ptm");
		putSourceCatParent("dehydroxylation", "ptm");
		putSourceCatParent("acetylation", "ptm");
		putSourceCatParent("deacetylation", "ptm");
		putSourceCatParent("glycosylation", "ptm");
		putSourceCatParent("deglycosylation", "ptm");
	}

	/**
	 * Provide some default pos/neg categories, like positive and negative regulation.
	 * Each of these categories should have at least one neutral (grand)parent!
	 */
	protected void insertDefaultContrasts()
	{
		addPosSourceCat("positive_regulation");
		addPosSourceCat("catalysis");
		addPosSourceCat("overexpression");

		addNegSourceCat("negative_regulation");
		addNegSourceCat("inhibition");
		addNegSourceCat("underexpression");

		addPosSourceCat("positive_genetic_interaction");

		addNegSourceCat("negative_genetic_interaction");
		addNegSourceCat("synthetic_lethality");
	}

	/**
	 * Provide a default mapping edge type to category mapping.
	 */
	protected void addDefaultPaintedParents()
	{
		/* SOURCE & CONSENSUS TYPES */
		
		sourceDraw.addColor("ptm", Color.BLUE, true);
		sourceDraw.addColor("colocalization", Color.YELLOW, true);
		sourceDraw.addColor("protein-dna_binding", Color.CYAN, true);
		sourceDraw.addColor("positive_regulation", Color.GREEN, true);
		sourceDraw.addColor("negative_regulation", Color.RED, true);
		sourceDraw.addColor("positive_genetic_interaction", Color.GREEN, true);
		sourceDraw.addColor("negative_genetic_interaction", Color.RED, true);

		sourceDraw.addArrowHead("regulation", ArrowHead.ARROW, true);
		sourceDraw.addArrowHead("positive_regulation", ArrowHead.ARROW, true);
		sourceDraw.addArrowHead("negative_regulation", ArrowHead.T, true);
		sourceDraw.addArrowHead("genetic_interaction", ArrowHead.NONE, true);

		sourceDraw.addArrowHead("ptm", ArrowHead.DIAMOND, true);
		sourceDraw.addArrowHead("protein-dna_binding", ArrowHead.ARROW, true);

		sourceDraw.addArrowHead("colocalization", ArrowHead.NONE, true);
		sourceDraw.addArrowHead("coexpression", ArrowHead.NONE, true);
		
		/* DIFFERENTIAL TYPES */
		
		diffDraw.addColor("increase_ppi", new Color(51,255,0), false);
		diffDraw.addColor("decrease_ppi", new Color(255,153,51), false);
		
		diffDraw.addColor("increases_regulation", new Color(0,153,0), false);
		diffDraw.addColor("decreases_regulation", new Color(204,0,0), false);
		
		diffDraw.addColor("increases_phosphorylation", new Color(0,153,255), false);
		diffDraw.addColor("decreases_phosphorylation", new Color(153,0,204), false);
		
		diffDraw.addColor("increases_dephosphorylation", new Color(255,0,204), false);
		diffDraw.addColor("decreases_dephosphorylation", new Color(0,51,255), false);
		
		diffDraw.addArrowHead("increase_ppi", ArrowHead.NONE, false);
		diffDraw.addArrowHead("decrease_ppi", ArrowHead.NONE, false);
		
		diffDraw.addArrowHead("increases_regulation", ArrowHead.ARROW, false);
		diffDraw.addArrowHead("decreases_regulation", ArrowHead.ARROW, false);
		
		diffDraw.addArrowHead("increases_phosphorylation", ArrowHead.DIAMOND, false);
		diffDraw.addArrowHead("decreases_phosphorylation", ArrowHead.DIAMOND, false);
		
		diffDraw.addArrowHead("increases_dephosphorylation", ArrowHead.DIAMOND, false);
		diffDraw.addArrowHead("decreases_dephosphorylation", ArrowHead.DIAMOND, false);
	}

	/**
	 * Provide a default mapping edge type to category mapping.
	 */
	protected void insertDefaultMappings()
	{
		boolean overwrite = false;

		// binding, coexpression, colocalization

		addSourceCategoryMapping("ppi", "ppi", overwrite);
		addSourceCategoryMapping("validated_ppi", "ppi", overwrite);
		addSourceCategoryMapping("protein-protein interaction", "ppi", overwrite);
		addSourceCategoryMapping("complex", "ppi", overwrite);
		addSourceCategoryMapping("complex formation", "ppi", overwrite);
		addSourceCategoryMapping("bind", "ppi", overwrite);
		addSourceCategoryMapping("binds", "ppi", overwrite);
		addSourceCategoryMapping("binding", "ppi", overwrite);

		addSourceCategoryMapping("protein-dna binding", "protein-dna_binding", overwrite);

		addSourceCategoryMapping("transcription", "transcription", overwrite);
		addSourceCategoryMapping("transcribes", "transcription", overwrite);
		addSourceCategoryMapping("transcribe", "transcription", overwrite);
		addSourceCategoryMapping("tf", "transcription", overwrite);
		addSourceCategoryMapping("transcription factor", "transcription", overwrite);
		addSourceCategoryMapping("transcription-factor binding", "transcription", overwrite);
		addSourceCategoryMapping("tf binding", "transcription", overwrite);

		addSourceCategoryMapping("coexpression", "coexpression", overwrite);
		addSourceCategoryMapping("coexpressed", "coexpression", overwrite);
		addSourceCategoryMapping("coexpresses", "coexpression", overwrite);
		addSourceCategoryMapping("coexpress", "coexpression", overwrite);
		addSourceCategoryMapping("coexpr", "coexpression", overwrite);
		addSourceCategoryMapping("coexpresses with", "coexpression", overwrite);

		addSourceCategoryMapping("colocalization", "colocalization", overwrite);
		addSourceCategoryMapping("colocalisation", "colocalization", overwrite);
		addSourceCategoryMapping("colocalized", "colocalization", overwrite);
		addSourceCategoryMapping("colocalised", "colocalization", overwrite);
		addSourceCategoryMapping("colocalizes", "colocalization", overwrite);
		addSourceCategoryMapping("colocalises", "colocalization", overwrite);
		addSourceCategoryMapping("colocalize", "colocalization", overwrite);
		addSourceCategoryMapping("colocalise", "colocalization", overwrite);
		addSourceCategoryMapping("colocalizing", "colocalization", overwrite);
		addSourceCategoryMapping("colocalising", "colocalization", overwrite);
		addSourceCategoryMapping("colocalizes with", "colocalization", overwrite);
		addSourceCategoryMapping("colocalises with", "colocalization", overwrite);
		addSourceCategoryMapping("coloc", "colocalization", overwrite);
		
		// PTM

		addSourceCategoryMapping("methylation", "methylation", overwrite);
		addSourceCategoryMapping("methylates", "methylation", overwrite);
		addSourceCategoryMapping("methylate", "methylation", overwrite);
		addSourceCategoryMapping("methylating", "methylation", overwrite);
		addSourceCategoryMapping("meth", "methylation", overwrite);

		addSourceCategoryMapping("demethylation", "demethylation", overwrite);
		addSourceCategoryMapping("demethylates", "demethylation", overwrite);
		addSourceCategoryMapping("demethylate", "demethylation", overwrite);
		addSourceCategoryMapping("demethylating", "demethylation", overwrite);
		addSourceCategoryMapping("demeth", "demethylation", overwrite);

		addSourceCategoryMapping("phosphorylation", "phosphorylation", overwrite);
		addSourceCategoryMapping("phosphorylates", "phosphorylation", overwrite);
		addSourceCategoryMapping("phosphorylate", "phosphorylation", overwrite);
		addSourceCategoryMapping("phosphorylating", "phosphorylation", overwrite);
		addSourceCategoryMapping("phosph", "phosphorylation", overwrite);
		addSourceCategoryMapping("phos", "phosphorylation", overwrite);

		addSourceCategoryMapping("dephosphorylation", "dephosphorylation", overwrite);
		addSourceCategoryMapping("dephosphorylates", "dephosphorylation", overwrite);
		addSourceCategoryMapping("dephosphorylate", "dephosphorylation", overwrite);
		addSourceCategoryMapping("dephosphorylating", "dephosphorylation", overwrite);
		addSourceCategoryMapping("dephosph", "dephosphorylation", overwrite);
		addSourceCategoryMapping("dephos", "dephosphorylation", overwrite);

		addSourceCategoryMapping("ubiquitination", "ubiquitination", overwrite);
		addSourceCategoryMapping("ubiquitinates", "ubiquitination", overwrite);
		addSourceCategoryMapping("ubiquitinate", "ubiquitination", overwrite);
		addSourceCategoryMapping("ubiquitinating", "ubiquitination", overwrite);
		addSourceCategoryMapping("ubiquitinylation", "ubiquitination", overwrite);
		addSourceCategoryMapping("ubiquitinylates", "ubiquitination", overwrite);
		addSourceCategoryMapping("ubiquitinylate", "ubiquitination", overwrite);
		addSourceCategoryMapping("ubiquitinylating", "ubiquitination", overwrite);

		addSourceCategoryMapping("deubiquitination", "deubiquitination", overwrite);
		addSourceCategoryMapping("deubiquitinates", "deubiquitination", overwrite);
		addSourceCategoryMapping("deubiquitinate", "deubiquitination", overwrite);
		addSourceCategoryMapping("deubiquitinating", "deubiquitination", overwrite);

		addSourceCategoryMapping("glycosylation", "glycosylation", overwrite);
		addSourceCategoryMapping("glycosylates", "glycosylation", overwrite);
		addSourceCategoryMapping("glycosylate", "glycosylation", overwrite);
		addSourceCategoryMapping("glycosylating", "glycosylation", overwrite);
		addSourceCategoryMapping("glyc", "glycosylation", overwrite);
		
		addSourceCategoryMapping("deglycosylation", "deglycosylation", overwrite);
		addSourceCategoryMapping("deglycosylates", "deglycosylation", overwrite);
		addSourceCategoryMapping("deglycosylate", "deglycosylation", overwrite);
		addSourceCategoryMapping("deglycosylating", "deglycosylation", overwrite);
		addSourceCategoryMapping("deglyc", "deglycosylation", overwrite);

		addSourceCategoryMapping("acetylation", "acetylation", overwrite);
		addSourceCategoryMapping("acetylates", "acetylation", overwrite);
		addSourceCategoryMapping("acetylate", "acetylation", overwrite);
		addSourceCategoryMapping("acetylating", "acetylation", overwrite);

		addSourceCategoryMapping("deacetylation", "deacetylation", overwrite);
		addSourceCategoryMapping("deacetylates", "deacetylation", overwrite);
		addSourceCategoryMapping("deacetylate", "deacetylation", overwrite);
		addSourceCategoryMapping("deacetylating", "deacetylation", overwrite);

		addSourceCategoryMapping("hydroxylation", "hydroxylation", overwrite);
		addSourceCategoryMapping("hydroxylates", "hydroxylation", overwrite);
		addSourceCategoryMapping("hydroxylate", "hydroxylation", overwrite);
		addSourceCategoryMapping("hydroxylating", "hydroxylation", overwrite);
		addSourceCategoryMapping("hydrox", "hydroxylation", overwrite);

		addSourceCategoryMapping("dehydroxylation", "dehydroxylation", overwrite);
		addSourceCategoryMapping("dehydroxylates", "dehydroxylation", overwrite);
		addSourceCategoryMapping("dehydroxylate", "dehydroxylation", overwrite);
		addSourceCategoryMapping("dehydroxylating", "dehydroxylation", overwrite);
		addSourceCategoryMapping("dehydrox", "dehydroxylation", overwrite);

		addSourceCategoryMapping("ptm", "ptm", overwrite);
		addSourceCategoryMapping("post-translational_modification", "ptm", overwrite);

		// regulation category and common synonyms
		addSourceCategoryMapping("regulation", "regulation", overwrite);
		addSourceCategoryMapping("regulates", "regulation", overwrite);
		addSourceCategoryMapping("regulate", "regulation", overwrite);
		addSourceCategoryMapping("regulating", "regulation", overwrite);
		addSourceCategoryMapping("influence", "regulation", overwrite);
		addSourceCategoryMapping("influences", "regulation", overwrite);
		addSourceCategoryMapping("influencing", "regulation", overwrite);
		addSourceCategoryMapping("effect", "regulation", overwrite);
		addSourceCategoryMapping("effects", "regulation", overwrite);
		addSourceCategoryMapping("effecting", "regulation", overwrite);
		addSourceCategoryMapping("affect", "regulation", overwrite);
		addSourceCategoryMapping("affects", "regulation", overwrite);
		addSourceCategoryMapping("affecting", "regulation", overwrite);
		addSourceCategoryMapping("unknown_regulation", "regulation", overwrite);
		addSourceCategoryMapping("gene regulation", "regulation", overwrite);
		addSourceCategoryMapping("protein regulation", "regulation", overwrite);

		// positive regulation category and common synonyms
		addSourceCategoryMapping("positive regulation", "positive_regulation", overwrite);
		addSourceCategoryMapping("positive reg", "positive_regulation", overwrite);
		addSourceCategoryMapping("pos regulation", "positive_regulation", overwrite);
		addSourceCategoryMapping("pos reg", "positive_regulation", overwrite);
		addSourceCategoryMapping("pos", "positive_regulation", overwrite);
		addSourceCategoryMapping("positive", "positive_regulation", overwrite);
		addSourceCategoryMapping("positively regulates", "positive_regulation", overwrite);
		addSourceCategoryMapping("positively regulating", "positive_regulation", overwrite);
		addSourceCategoryMapping("upregulation", "positive_regulation", overwrite);
		addSourceCategoryMapping("upregulated", "positive_regulation", overwrite);
		addSourceCategoryMapping("activate", "positive_regulation", overwrite);
		addSourceCategoryMapping("activates", "positive_regulation", overwrite);
		addSourceCategoryMapping("activating", "positive_regulation", overwrite);
		addSourceCategoryMapping("activation", "positive_regulation", overwrite);

		addSourceCategoryMapping("catalysis", "catalysis", overwrite);
		addSourceCategoryMapping("catalyzis", "catalysis", overwrite);
		addSourceCategoryMapping("catalyses", "catalysis", overwrite);
		addSourceCategoryMapping("catalyzes", "catalysis", overwrite);
		addSourceCategoryMapping("catalysing", "catalysis", overwrite);
		addSourceCategoryMapping("catalyzing", "catalysis", overwrite);
		addSourceCategoryMapping("catalysation", "catalysis", overwrite);
		addSourceCategoryMapping("catalyzation", "catalysis", overwrite);

		// negative regulation category and common synonyms
		addSourceCategoryMapping("negative regulation", "negative_regulation", overwrite);
		addSourceCategoryMapping("negative reg", "negative_regulation", overwrite);
		addSourceCategoryMapping("neg regulation", "negative_regulation", overwrite);
		addSourceCategoryMapping("neg reg", "negative_regulation", overwrite);
		addSourceCategoryMapping("neg", "negative_regulation", overwrite);
		addSourceCategoryMapping("negative", "negative_regulation", overwrite);
		addSourceCategoryMapping("negatively regulates", "negative_regulation", overwrite);
		addSourceCategoryMapping("negatively regulating", "negative_regulation", overwrite);
		addSourceCategoryMapping("downregulation", "negative_regulation", overwrite);
		addSourceCategoryMapping("downregulated", "negative_regulation", overwrite);
		addSourceCategoryMapping("inactivation", "negative_regulation", overwrite);
		addSourceCategoryMapping("inactivate", "negative_regulation", overwrite);
		addSourceCategoryMapping("inactivates", "negative_regulation", overwrite);
		addSourceCategoryMapping("inactivated", "negative_regulation", overwrite);

		addSourceCategoryMapping("inhibit", "inhibition", overwrite);
		addSourceCategoryMapping("inhibits", "inhibition", overwrite);
		addSourceCategoryMapping("inhibiting", "inhibition", overwrite);
		addSourceCategoryMapping("inhibition", "inhibition", overwrite);
		addSourceCategoryMapping("repress", "inhibition", overwrite);
		addSourceCategoryMapping("represses", "inhibition", overwrite);
		addSourceCategoryMapping("repressing", "inhibition", overwrite);
		addSourceCategoryMapping("repression", "inhibition", overwrite);
		
		//  (virtual) expression edges
		
		addSourceCategoryMapping("expression", "expression", overwrite);
		addSourceCategoryMapping("expresses", "expression", overwrite);
		addSourceCategoryMapping("expressed", "expression", overwrite);
		addSourceCategoryMapping("express", "expression", overwrite);
		
		addSourceCategoryMapping("overexpression", "overexpression", overwrite);
		addSourceCategoryMapping("overexpresses", "overexpression", overwrite);
		addSourceCategoryMapping("overexpressed", "overexpression", overwrite);
		addSourceCategoryMapping("overexpress", "overexpression", overwrite);
		
		addSourceCategoryMapping("underexpression", "underexpression", overwrite);
		addSourceCategoryMapping("underexpresses", "underexpression", overwrite);
		addSourceCategoryMapping("underexpressed", "underexpression", overwrite);
		addSourceCategoryMapping("underexpress", "underexpression", overwrite);

		// genetic interactions
		addSourceCategoryMapping("positive genetic interaction", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("positive gi", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("alleviating gi", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("alleviating genetic interaction", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("positive epistatic interaction", "positive_genetic_interaction", overwrite);
		addSourceCategoryMapping("alleviating epistatic interaction", "positive_genetic_interaction", overwrite);

		addSourceCategoryMapping("negative genetic interaction", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("negative gi", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("aggrevating genetic interaction", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("aggrevating gi", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("negative epistatic interaction", "negative_genetic_interaction", overwrite);
		addSourceCategoryMapping("aggrevating epistatic interaction", "negative_genetic_interaction", overwrite);

		addSourceCategoryMapping("synthetic lethality", "synthetic_lethality", overwrite);
		addSourceCategoryMapping("synth lethality", "synthetic_lethality", overwrite);
		addSourceCategoryMapping("synthetically lethal", "synthetic_lethality", overwrite);
		addSourceCategoryMapping("synth lethal", "synthetic_lethality", overwrite);
		addSourceCategoryMapping("sl", "synthetic_lethality", overwrite);

		addSourceCategoryMapping("genetic_interaction", "genetic_interaction", overwrite);
		addSourceCategoryMapping("genetically_interacting", "genetic_interaction", overwrite);
		addSourceCategoryMapping("gi", "genetic_interaction", overwrite);
		addSourceCategoryMapping("epistatic", "genetic_interaction", overwrite);
		addSourceCategoryMapping("epistasis", "genetic_interaction", overwrite);
		addSourceCategoryMapping("epistatic interaction", "genetic_interaction", overwrite);
	}

}
