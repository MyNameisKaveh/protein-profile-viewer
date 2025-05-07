# Imports (requests, gradio, Bio, Counter, matplotlib, io, PIL, sys, traceback, json)
import gradio as gr
import requests
from Bio.SeqUtils import molecular_weight
from collections import Counter
import matplotlib.pyplot as plt
import io
from PIL import Image
import sys 
import traceback
import json 

# Constants
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
AMINO_ACID_NAMES = {
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic acid',
    'C': 'Cysteine', 'Q': 'Glutamine', 'E': 'Glutamic acid', 'G': 'Glycine',
    'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
    'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
    'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine'
}
STANDARD_AMINO_ACIDS_ORDER = "ARNDCQEGHILKMFPSTWYV" 

# --- Helper Functions (Extraction & Plotting) ---

def get_amino_acid_frequencies(sequence):
    # ... (Implementation from previous version) ...
    if not sequence or sequence == "N/A": return None, "Sequence not available for analysis."
    cleaned_sequence = "".join(filter(lambda x: x in AMINO_ACID_NAMES, sequence.upper()))
    if not cleaned_sequence: return None, "No valid amino acids found in sequence for counting."
    counts = Counter(cleaned_sequence)
    frequencies = {aa: counts.get(aa, 0) for aa in STANDARD_AMINO_ACIDS_ORDER}
    return frequencies, None

def plot_amino_acid_frequencies(frequencies):
    # ... (Implementation from previous version) ...
    if not frequencies: return None
    ordered_keys = [key for key in STANDARD_AMINO_ACIDS_ORDER if key in frequencies]
    labels = [f"{aa}: {AMINO_ACID_NAMES.get(aa, aa)}" for aa in ordered_keys]
    values = [frequencies[aa] for aa in ordered_keys]
    fig, ax = plt.subplots(figsize=(12, 7)); ax.bar(labels, values, color='skyblue')
    ax.set_xlabel("Amino Acid"); ax.set_ylabel("Frequency"); ax.set_title("Amino Acid Frequency Plot")
    plt.xticks(rotation=75, ha="right", fontsize=8); plt.tight_layout()
    buf = io.BytesIO(); plt.savefig(buf, format='png'); buf.seek(0)
    img = Image.open(buf); plt.close(fig)
    return img

def extract_sequence_features(uniprot_data):
    # ... (Implementation from previous version) ...
    features_of_interest_uppercase = {
        "DOMAIN": "blue", "MOTIF": "green", "ACTIVE_SITE": "red",
        "BINDING_SITE": "orange", "MOD_RES": "purple", 
        "HELIX": "cyan", "STRAND": "magenta", "TURN": "gold"
    }
    extracted_features = []
    if "features" in uniprot_data and uniprot_data["features"] is not None:
        for feature_item in uniprot_data["features"]:
            feature_type_raw = feature_item.get("type")
            if not isinstance(feature_type_raw, str): continue
            feature_type_normalized = feature_type_raw.strip().upper()
            if feature_type_normalized in features_of_interest_uppercase:
                try:
                    location_obj = feature_item.get("location", {})
                    begin_pos_val_str = None; end_pos_val_str = None
                    start_node = location_obj.get("start"); end_node = location_obj.get("end")
                    position_node = location_obj.get("position")
                    if start_node and isinstance(start_node, dict) and "value" in start_node: begin_pos_val_str = str(start_node["value"])
                    if end_node and isinstance(end_node, dict) and "value" in end_node: end_pos_val_str = str(end_node["value"])
                    if position_node and isinstance(position_node, dict) and "value" in position_node:
                        pos_val_str = str(position_node["value"])
                        if begin_pos_val_str is None: begin_pos_val_str = pos_val_str
                        if end_pos_val_str is None: end_pos_val_str = pos_val_str
                        if "start" not in location_obj and "end" not in location_obj:
                             begin_pos_val_str = pos_val_str; end_pos_val_str = pos_val_str
                    if begin_pos_val_str is None or end_pos_val_str is None: continue
                    begin_pos = int(begin_pos_val_str); end_pos = int(end_pos_val_str)
                    if begin_pos > end_pos: continue
                    extracted_features.append({
                        "type": feature_type_raw, "begin": begin_pos, "end": end_pos,
                        "description": feature_item.get("description", feature_type_raw), 
                        "color": features_of_interest_uppercase[feature_type_normalized]
                    })
                except (ValueError, TypeError, AttributeError): continue 
    return extracted_features

def plot_sequence_features(sequence_length, features):
    # ... (Implementation from previous version) ...
    if not features or sequence_length == 0: return None
    fig, ax = plt.subplots(figsize=(12, max(3, len(features) * 0.4) + 1.5)) 
    ax.set_xlim(0, sequence_length); ax.set_xlabel("Amino Acid Position"); ax.set_yticks([]) 
    ax.set_title("Sequence Features Plot"); legend_handles = {} 
    y_pos_counter = 0; plotted_feature_types_in_legend = set(); bar_height = 0.8 
    for feature in sorted(features, key=lambda x: x["begin"]):
        begin = feature["begin"]; end = feature["end"]; color = feature["color"]
        width = max(1, end - begin + 1) 
        ax.barh(y_pos_counter, width, height=bar_height, left=begin -1, color=color, edgecolor='black', alpha=0.7)
        if feature['type'] not in plotted_feature_types_in_legend:
            legend_handles[feature['type']] = plt.Rectangle((0, 0), 1, 1, fc=color, alpha=0.7)
            plotted_feature_types_in_legend.add(feature['type'])
        y_pos_counter += 1 
    if y_pos_counter > 0: ax.set_ylim(-0.5, y_pos_counter -1 + bar_height/2 + 0.5) 
    else: plt.close(fig); return None
    if legend_handles: ax.legend(legend_handles.values(), legend_handles.keys(), title="Feature Types", bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    plt.tight_layout(rect=[0, 0, 0.83, 0.96])
    buf = io.BytesIO(); plt.savefig(buf, format='png'); buf.seek(0)
    img = Image.open(buf); plt.close(fig)
    return img

def extract_interactions(uniprot_data):
    # ... (Implementation from previous version) ...
    interactions_list = []
    if "comments" in uniprot_data:
        for comment in uniprot_data["comments"]:
            if comment.get("commentType") == "INTERACTION" and "interactions" in comment:
                for interaction_entry in comment["interactions"]:
                    interactant_one_acc = interaction_entry.get("interactantOne", {}).get("uniProtKBAccession")
                    interactant_two_acc = interaction_entry.get("interactantTwo", {}).get("uniProtKBAccession")
                    interactant_two_gene = interaction_entry.get("interactantTwo", {}).get("geneName")
                    if interactant_one_acc == uniprot_data.get("primaryAccession") and interactant_two_acc:
                        partner_display_name = interactant_two_acc
                        if interactant_two_gene:
                            partner_display_name = f"{interactant_two_gene} ({interactant_two_acc})"
                        interactions_list.append(f"- Interacts with: **{partner_display_name}**")
    if not interactions_list: return "No specific interaction partners listed in UniProt comments."
    return "\n".join(sorted(interactions_list))

def extract_pathways(uniprot_data):
    # ... (Implementation from previous version) ...
    pathways = []
    pathway_databases = { "KEGG": "https://www.genome.jp/dbget-bin/www_bget?", "Reactome": "https://reactome.org/content/detail/" }
    if "uniProtKBCrossReferences" in uniprot_data:
        for xref in uniprot_data["uniProtKBCrossReferences"]:
            db_name = xref.get("database")
            if db_name in pathway_databases:
                pathway_id = xref.get("id"); pathway_description = ""
                if "properties" in xref:
                    for prop in xref["properties"]:
                        if prop.get("key") == "PathwayName" or prop.get("key") == "Description":
                            pathway_description = prop.get("value"); break
                if pathway_id:
                    link = pathway_databases[db_name] + pathway_id
                    display_text = f"{pathway_description} ({pathway_id})" if pathway_description else pathway_id
                    pathways.append(f"- [{display_text}]({link}) ({db_name})")
    if not pathways: return "No pathway information found in KEGG or Reactome cross-references."
    return "\n".join(sorted(list(set(pathways))))

def extract_disease_info(uniprot_data):
    # ... (Implementation from previous version) ...
    disease_info_list = []
    if "comments" in uniprot_data:
        for comment in uniprot_data["comments"]:
            if comment.get("commentType") == "DISEASE" and "disease" in comment:
                disease_entry = comment["disease"]
                disease_name = disease_entry.get("diseaseId", "Unknown disease")
                description = disease_entry.get("description", "No description available.")
                mim_id = None
                if "diseaseCrossReference" in disease_entry and disease_entry["diseaseCrossReference"].get("database") == "MIM":
                    mim_id = disease_entry["diseaseCrossReference"].get("id")
                disease_md = f"**{disease_name}**"
                if mim_id: disease_md += f" (MIM: [{mim_id}](https://www.omim.org/entry/{mim_id}))"
                disease_md += f"\n   - *Description:* {description}\n"
                if "note" in comment and "texts" in comment["note"]:
                    for note_text_obj in comment["note"]["texts"]:
                        note_val = note_text_obj.get("value")
                        if note_val: disease_md += f"   - *Note:* {note_val}\n"
                disease_info_list.append(disease_md)
    if not disease_info_list: return "No specific disease association information found in UniProt comments."
    return "\n---\n".join(sorted(disease_info_list))

def extract_publications(uniprot_data):
    # ... (Implementation from previous version) ...
    publications_list = []
    if "references" in uniprot_data:
        for ref_idx, ref in enumerate(uniprot_data.get("references", [])):
            citation = ref.get("citation", {})
            title = citation.get("title", "N/A")
            authors = ", ".join(citation.get("authors", ["N/A"]))
            journal = citation.get("journalName", "N/A"); volume = citation.get("volume", "")
            first_page = citation.get("firstPage", ""); last_page = citation.get("lastPage", "")
            publication_date = citation.get("publicationDate", "")
            pubmed_id = None; doi_id = None
            if "citationCrossReferences" in citation:
                for xref_cite in citation["citationCrossReferences"]:
                    if xref_cite.get("database") == "PubMed": pubmed_id = xref_cite.get("id")
                    elif xref_cite.get("database") == "DOI": doi_id = xref_cite.get("id")
            pub_md = f"**{ref_idx + 1}. Title:** {title}\n"
            pub_md += f"   - *Authors:* {authors}\n   - *Journal:* {journal}"
            if volume: pub_md += f", Vol. {volume}"
            if first_page: pub_md += f", pp. {first_page}"
            if last_page: pub_md += f"-{last_page}"
            if publication_date: pub_md += f" ({publication_date})\n" else: pub_md += "\n"
            if pubmed_id: pub_md += f"   - *PubMed:* [{pubmed_id}](https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/)\n"
            if doi_id: pub_md += f"   - *DOI:* [{doi_id}](https://doi.org/{doi_id})\n"
            publications_list.append(pub_md)
    if not publications_list: return "No publication information found in this UniProt entry."
    return "\n---\n".join(publications_list)

def extract_cross_references(uniprot_data):
    # ... (Implementation from previous version) ...
    xref_list = []
    target_databases = {
        "Ensembl": "https://www.ensembl.org/id/", "GeneID": "https://www.ncbi.nlm.nih.gov/gene/", 
        "RefSeq": "https://www.ncbi.nlm.nih.gov/nuccore/", "GO": "https://amigo.geneontology.org/amigo/term/", 
        "InterPro": "https://www.ebi.ac.uk/interpro/entry/InterPro/", 
        "Pfam": "https://www.ebi.ac.uk/interpro/entry/pfam/", "PDB": "https://www.rcsb.org/structure/",
        "KEGG": "https://www.genome.jp/dbget-bin/www_bget?", "Reactome": "https://reactome.org/content/detail/"
    }
    grouped_xrefs = {db: [] for db in target_databases}
    if "uniProtKBCrossReferences" in uniprot_data:
        for xref in uniprot_data["uniProtKBCrossReferences"]:
            db_name = xref.get("database")
            if db_name in target_databases:
                xref_id = xref.get("id"); link_url = None; display_text = xref_id
                if db_name == "Ensembl" and "properties" in xref:
                    for prop in xref["properties"]:
                        if prop.get("key") == "GeneId":
                            ensembl_gene_id = prop.get("value")
                            link_url = f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={ensembl_gene_id}"; display_text = ensembl_gene_id; break
                    if link_url is None:
                        for prop in xref["properties"]:
                           if prop.get("key") == "ProteinId":
                                ensembl_prot_id = prop.get("value")
                                link_url = f"https://www.ensembl.org/Homo_sapiens/Transcript/Summary?p={ensembl_prot_id}"; display_text = ensembl_prot_id; break
                elif db_name == "RefSeq" and "properties" in xref:
                     for prop in xref["properties"]:
                         if prop.get("key") == "ProteinId" or prop.get("key") == "NucleotideSequenceId":
                             refseq_id = prop.get("value"); link_url = target_databases[db_name] + refseq_id; display_text = refseq_id; break
                elif db_name == "GO" and xref_id:
                    term_name = xref_id
                    if "properties" in xref:
                        for prop in xref["properties"]:
                             if prop.get("key") == "GoTerm": term_name = prop.get("value"); break
                    link_url = target_databases[db_name] + xref_id; display_text = f"{term_name} ({xref_id})"
                elif db_name == "Pfam" and xref_id: link_url = f"https://www.ebi.ac.uk/interpro/entry/pfam/{xref_id}"
                elif xref_id: link_url = target_databases[db_name] + xref_id
                if link_url: grouped_xrefs[db_name].append(f"[{display_text}]({link_url})")
        for db_name, links in grouped_xrefs.items():
            if links: xref_list.append(f"**{db_name}:** " + ", ".join(sorted(list(set(links)))))
    if not xref_list: return "No cross-references found for the selected databases."
    return "\n".join(xref_list)

# --- NEW Search Function ---
def search_uniprot_by_name(search_term, result_limit=5):
    if not search_term or len(search_term) < 3: return "Enter at least 3 characters."
    params = { "query": f'({search_term}) AND (reviewed:true)', "fields": "accession,protein_name,organism_name,id", "format": "json", "size": result_limit }
    md_parts = [f"### Search Results for '{search_term}' (Top {result_limit} Reviewed):\n"]
    try:
        response = requests.get(UNIPROT_SEARCH_URL, params=params); response.raise_for_status(); data = response.json()
        results = data.get("results")
        if not results: md_parts.append("No reviewed entries found.")
        else:
            md_parts.append(f"*Found {len(results)}. Copy desired Accession ID and paste below.*\n---")
            for entry in results:
                acc = entry.get("primaryAccession", "N/A"); protein_name = "N/A"
                if entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"):
                    protein_name = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
                elif entry.get("proteinDescription", {}).get("submissionNames"): protein_name = entry["proteinDescription"]["submissionNames"][0].get("fullName",{}).get("value", "N/A")
                organism = entry.get("organism", {}).get("scientificName", "N/A"); uniprot_id_display = entry.get("uniProtkbId", "")
                md_parts.append(f"- **{acc}** (`{uniprot_id_display}`): {protein_name} - *{organism}*")
            md_parts.append("---")
    except Exception as e: md_parts.append(f"\n**Search Error:** {e}")
    return "\n".join(md_parts)

# --- Main Function to Get All Info ---
def get_protein_info(uniprot_id):
    empty_plot = gr.update(value=None, visible=False); empty_str = ""; err_msg = "Error"
    outputs_on_error = (err_msg, empty_plot, empty_plot, err_msg, err_msg, err_msg, err_msg, err_msg, err_msg) 
    if not uniprot_id: return ("Enter UniProt ID.",) + outputs_on_error[1:]
    
    url = UNIPROT_API_URL.format(accession=uniprot_id.strip().upper())
    try:
        response = requests.get(url); response.raise_for_status(); data = response.json()
        acc = data.get("primaryAccession", "N/A"); id_display = data.get("uniProtkbId", "N/A")
        link = f"https://www.uniprot.org/uniprotkb/{acc}/entry"
        name_dict = data.get("proteinDescription", {}).get("recommendedName", {}); name = name_dict.get("fullName", {}).get("value", "N/A")
        if name == "N/A" and data.get("proteinDescription", {}).get("submissionNames"): name = data["proteinDescription"]["submissionNames"][0].get("fullName", {}).get("value", "N/A")
        genes_data = data.get("genes"); gene_str = "N/A"
        if genes_data:
            g_list = [g.get("geneName", {}).get("value", "") for g in genes_data if g.get("geneName")]
            if not g_list: g_list = [g.get("orfNames", [{}])[0].get("value", "") for g in genes_data if g.get("orfNames")]
            gene_str = ", ".join(filter(None, g_list)) or "N/A"
        org_data = data.get("organism", {}); org_sci_name = org_data.get("scientificName", "N/A"); org_common_name = org_data.get("commonName", "")
        org_display = f"{org_sci_name}" + (f" ({org_common_name})" if org_common_name else "")
        seq_info = data.get("sequence", {}); seq_val = seq_info.get("value", "N/A"); length = seq_info.get("length", 0)
        status_val = data.get("entryAudit", {}).get("entryType", "N/A").replace("UniProtKB ", "")
        existence_val = data.get("proteinExistence", "N/A").replace(": Evidence at ", ": ")
        score_val = data.get("annotationScore", "N/A"); mw_str = "N/A"
        if seq_val != "N/A" and length > 0:
            try:
                clean_seq = "".join(filter(lambda x: x in AMINO_ACID_NAMES, seq_val.upper()))
                if clean_seq: mw_str = f"{molecular_weight(clean_seq, seq_type='protein'):.2f} Da"
                else: mw_str = "Invalid sequence for MW"
            except: mw_str = "Error in MW calc"
        comments_data = data.get("comments", []); func_comment = "N/A" 
        for c_item in comments_data:
            if c_item.get("commentType") == "FUNCTION":
                texts = c_item.get("texts", []);
                if texts: func_comment = texts[0].get("value", "N/A"); break
        
        overview_md = (f"## {id_display} ({acc})\n[{acc} on UniProt]({link})\n\n"
                       f"**Protein:** {name}\n**Gene:** {gene_str}\n**Status:** {status_val}\n"
                       f"**Organism:** {org_display}\n**Length:** {length} aa\n"
                       f"**Existence:** {existence_val}\n**Score:** {score_val}/5\n**Calc. MW:** {mw_str}\n\n"
                       f"**Function Snippet:**\n{func_comment}\n\n"
                       f"**Sequence (first 100 aa):**\n`{seq_val[:100]}{'...' if len(seq_val) > 100 else ''}`\n\n"
                       f"--- \n*More details in other tabs.*")
        
        interactions_md = extract_interactions(data); pathways_md = extract_pathways(data)
        disease_md = extract_disease_info(data); publications_md = extract_publications(data)
        xref_md = extract_cross_references(data) 

        aa_freq, aa_err = get_amino_acid_frequencies(seq_val); aa_plot_upd = empty_plot
        if aa_err: overview_md += f"\n\n**AA Freq Error:** {aa_err}"
        elif aa_freq:
            img_aa = plot_amino_acid_frequencies(aa_freq)
            if img_aa: aa_plot_upd = gr.update(value=img_aa, visible=True)
        
        seq_feat = extract_sequence_features(data); feat_plot_upd = empty_plot; feat_msg = ""
        if seq_feat and length > 0:
            img_feat = plot_sequence_features(length, seq_feat)
            if img_feat: feat_plot_upd = gr.update(value=img_feat, visible=True)
            else: feat_msg = "Could not generate feature plot."
        elif not seq_feat and length > 0 : feat_msg = "No relevant features found for plotting."
        
        return (overview_md, aa_plot_upd, feat_plot_upd, feat_msg, 
                pathways_md, interactions_md, disease_md, publications_md, xref_md)

    except requests.exceptions.HTTPError as e: 
        err_msg_http = f"Error: ID '{uniprot_id}' not found." if e.response.status_code == 404 else f"HTTP error: {e}"
        return (err_msg_http,) + outputs_on_error[1:]
    except Exception as e:
        return (f"Error: {str(e)[:150]}",) + outputs_on_error[1:]

# --- Gradio UI Definition ---
with gr.Blocks(theme=gr.themes.Glass()) as iface:
    gr.Markdown("# Protein Profile Viewer (v1.5 - ID Search)")
    gr.Markdown("Enter a UniProt ID directly, **OR** search for a protein/gene name below to find its ID first.")
    
    with gr.Group():
        gr.Markdown("### Find UniProt ID by Name/Keyword")
        with gr.Row():
            search_term_input = gr.Textbox(label="Enter Protein/Gene Name (e.g., insulin, EGFR)", scale=3)
            search_button = gr.Button("Search ID", scale=1)
        search_results_output = gr.Markdown(label="Search Results (Copy the Accession ID you need)")

    gr.Markdown("---") 

    gr.Markdown("### View Protein Profile")
    with gr.Row():
        protein_id_input = gr.Textbox(label="Enter UniProt Accession ID (e.g., P01308 or P00533)", placeholder="Paste ID here", scale=3)
        submit_button = gr.Button("Get Profile", scale=1, variant="primary")

    with gr.Tabs():
        with gr.TabItem("Overview"):
            gr.Markdown("### Protein Overview\nKey info: UniProt ID, name, gene, organism, function, sequence snippet & UniProt link.")
            overview_output = gr.Markdown()
        with gr.TabItem("Analysis Plots"):
            gr.Markdown("### Sequence Analysis Visualizations\nAmino acid composition and annotated features (domains, motifs, etc.).")
            with gr.Column(): 
                aa_freq_plot_output = gr.Image(label="Amino Acid Frequency Plot", type="pil", show_label=True, visible=False)
                seq_features_plot_output = gr.Image(label="Sequence Features Plot", type="pil", show_label=True, visible=False)
                seq_features_message_output = gr.Markdown() 
        with gr.TabItem("Functional Context"):
            gr.Markdown("### Pathways, Interactions & Disease\nBiological context: pathways, interaction partners, and associated diseases.")
            with gr.Accordion("Biological Pathways (KEGG, Reactome)", open=False):
                 pathways_output = gr.Markdown()
            with gr.Accordion("Protein Interactions", open=False):
                 interactions_output = gr.Markdown()
            with gr.Accordion("Disease Associations", open=False):
                 disease_output = gr.Markdown()
        with gr.TabItem("Publications"):
            gr.Markdown("### Relevant Publications\nA list of scientific publications from UniProt.")
            publications_output = gr.Markdown()
        with gr.TabItem("Cross-references"):
             gr.Markdown("### Database Links\nLinks to this protein's entry in other relevant biological databases.")
             xref_output = gr.Markdown()

    search_button.click(fn=search_uniprot_by_name, inputs=search_term_input, outputs=search_results_output)
    submit_button.click(
        fn=get_protein_info,
        inputs=protein_id_input,
        outputs=[overview_output, aa_freq_plot_output, seq_features_plot_output, seq_features_message_output, 
                 pathways_output, interactions_output, disease_output, publications_output, xref_output] 
    )
    gr.Examples(
        examples=[["P05067"], ["P00533"], ["Q9BYF1"], ["P0DP23"], ["P04637"]],
        inputs=protein_id_input 
    )

if __name__ == "__main__":
    iface.launch()
