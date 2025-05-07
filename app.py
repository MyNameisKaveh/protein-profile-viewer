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

UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"

AMINO_ACID_NAMES = {
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic acid',
    'C': 'Cysteine', 'Q': 'Glutamine', 'E': 'Glutamic acid', 'G': 'Glycine',
    'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
    'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
    'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine'
}
STANDARD_AMINO_ACIDS_ORDER = "ARNDCQEGHILKMFPSTWYV" 

def get_amino_acid_frequencies(sequence):
    if not sequence or sequence == "N/A": return None, "Sequence not available for analysis."
    cleaned_sequence = "".join(filter(lambda x: x in AMINO_ACID_NAMES, sequence.upper()))
    if not cleaned_sequence: return None, "No valid amino acids found in sequence for counting."
    counts = Counter(cleaned_sequence)
    frequencies = {aa: counts.get(aa, 0) for aa in STANDARD_AMINO_ACIDS_ORDER}
    return frequencies, None

def plot_amino_acid_frequencies(frequencies):
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
    if not interactions_list:
        return "No specific interaction partners listed in UniProt comments."
    return "\n".join(interactions_list)

def extract_pathways(uniprot_data):
    pathways = []
    pathway_databases = {
        "KEGG": "https://www.genome.jp/dbget-bin/www_bget?",
        "Reactome": "https://reactome.org/content/detail/"
    }
    if "uniProtKBCrossReferences" in uniprot_data:
        for xref in uniprot_data["uniProtKBCrossReferences"]:
            db_name = xref.get("database")
            if db_name in pathway_databases:
                pathway_id = xref.get("id")
                pathway_description = ""
                if "properties" in xref:
                    for prop in xref["properties"]:
                        if prop.get("key") == "PathwayName" or prop.get("key") == "Description":
                            pathway_description = prop.get("value")
                            break
                if pathway_id:
                    link = pathway_databases[db_name] + pathway_id
                    display_text = f"{pathway_description} ({pathway_id})" if pathway_description else pathway_id
                    pathways.append(f"- [{display_text}]({link}) ({db_name})")
    if not pathways:
        return "No pathway information found in KEGG or Reactome cross-references."
    return "\n".join(sorted(list(set(pathways))))

def get_protein_info(uniprot_id):
    empty_plot_update = gr.update(value=None, visible=False); empty_str = ""
    default_error_msg = "Error fetching or processing data."

    if not uniprot_id:
        return ("Please enter a UniProt ID.", empty_plot_update, empty_plot_update, 
                empty_str, "No data to display.", "No data to display.")
    
    url = UNIPROT_API_URL.format(accession=uniprot_id.strip().upper())
    try:
        response = requests.get(url); response.raise_for_status(); uniprot_api_data = response.json()
        
        primary_accession = uniprot_api_data.get("primaryAccession", "N/A")
        uniprot_link = f"https://www.uniprot.org/uniprotkb/{primary_accession}/entry"
        
        protein_data_dict = uniprot_api_data.get("proteinDescription", {}).get("recommendedName", {})
        protein_name = protein_data_dict.get("fullName", {}).get("value", "N/A")
        if protein_name == "N/A" and uniprot_api_data.get("proteinDescription", {}).get("submissionNames"):
            protein_name = uniprot_api_data["proteinDescription"]["submissionNames"][0].get("fullName", {}).get("value", "N/A")
        gene_names_data = uniprot_api_data.get("genes"); gene_name_str = "N/A"
        if gene_names_data:
            gene_names_list = [g.get("geneName", {}).get("value", "") for g in gene_names_data if g.get("geneName")]
            if not gene_names_list: gene_names_list = [g.get("orfNames", [{}])[0].get("value", "") for g in gene_names_data if g.get("orfNames")]
            gene_name_str = ", ".join(filter(None, gene_names_list)) or "N/A"
        organism_data = uniprot_api_data.get("organism", {}); organism_name = organism_data.get("scientificName", "N/A")
        if organism_data.get("commonName"): organism_name += f" ({organism_data.get('commonName')})"
        sequence_info = uniprot_api_data.get("sequence", {}); sequence = sequence_info.get("value", "N/A"); length = sequence_info.get("length", 0)
        calculated_mol_weight_str = "N/A"
        if sequence != "N/A" and length > 0:
            try:
                cleaned_sequence_for_mw = "".join(filter(lambda x: x in AMINO_ACID_NAMES, sequence.upper()))
                if cleaned_sequence_for_mw:
                     calculated_mol_weight = molecular_weight(cleaned_sequence_for_mw, seq_type="protein")
                     calculated_mol_weight_str = f"{calculated_mol_weight:.2f} Da"
                else: calculated_mol_weight_str = "Invalid sequence for MW calculation"
            except Exception: calculated_mol_weight_str = "Error in MW calculation"
        comments = uniprot_api_data.get("comments", []); function_comment_en = "N/A" 
        for comment_item in comments:
            if comment_item.get("commentType") == "FUNCTION":
                function_texts = comment_item.get("texts", [])
                if function_texts: function_comment_en = function_texts[0].get("value", "N/A"); break
        
        overview_md = (
            f"**Primary Accession:** {primary_accession} ([View on UniProt]({uniprot_link}))\n" # Added UniProt link
            f"**Protein Name:** {protein_name}\n"
            f"**Gene Name(s):** {gene_name_str}\n"
            f"**Organism:** {organism_name}\n"
            f"**Sequence Length:** {length} amino acids\n"
            f"**Molecular Weight (calculated):** {calculated_mol_weight_str}\n\n"
            f"**Function Snippet:**\n{function_comment_en}\n\n"
            f"**Amino Acid Sequence (first 100 residues):**\n`{sequence[:100]}{'...' if len(sequence) > 100 else ''}`"
        )
        
        interactions_md = extract_interactions(uniprot_api_data)
        pathways_md = extract_pathways(uniprot_api_data)

        aa_frequencies, freq_error = get_amino_acid_frequencies(sequence)
        aa_plot_update = empty_plot_update
        if freq_error: overview_md += f"\n\n**AA Freq. Analysis Error:** {freq_error}"
        elif aa_frequencies:
            pil_image_aa = plot_amino_acid_frequencies(aa_frequencies)
            if pil_image_aa: aa_plot_update = gr.update(value=pil_image_aa, visible=True)
        
        seq_features = extract_sequence_features(uniprot_api_data) 
        feature_plot_update = empty_plot_update; features_plot_message = ""
        if seq_features and length > 0:
            pil_image_features = plot_sequence_features(length, seq_features)
            if pil_image_features: feature_plot_update = gr.update(value=pil_image_features, visible=True)
            else: features_plot_message = "Could not generate sequence feature plot."
        elif not seq_features and length > 0 : features_plot_message = "No relevant features found for plotting."
        
        return overview_md, aa_plot_update, feature_plot_update, features_plot_message, pathways_md, interactions_md

    except requests.exceptions.HTTPError as http_err: 
        status_code = http_err.response.status_code if http_err.response is not None else "N/A"
        err_msg = f"Error: Protein ID '{uniprot_id}' not found." if status_code == 404 else f"HTTP error: {http_err}"
        return err_msg, empty_plot_update, empty_plot_update, default_error_msg, default_error_msg, default_error_msg
    except requests.exceptions.RequestException as req_err:
        return f"Network error: {req_err}", empty_plot_update, empty_plot_update, default_error_msg, default_error_msg, default_error_msg
    except Exception as e:
        return f"Unexpected error: {str(e)[:150]}", empty_plot_update, empty_plot_update, default_error_msg, default_error_msg, default_error_msg

with gr.Blocks(theme=gr.themes.Soft()) as iface:
    gr.Markdown("# Protein Profile Viewer (v1.1 - UI/Data Order Impr.)")
    gr.Markdown("Enter a UniProt ID to display its basic information, sequence analysis, pathways, and interactions.")
    
    with gr.Row():
        protein_id_input = gr.Textbox(label="Enter UniProt ID (e.g., P05067 or P00533)", placeholder="e.g., P0DP23", scale=3)
        submit_button = gr.Button("Submit", scale=1, variant="primary")

    with gr.Tabs():
        with gr.TabItem("Overview & Sequence"):
            overview_output = gr.Markdown(label="Protein Information")
        with gr.TabItem("Sequence Analysis Plots"):
            with gr.Column(): 
                aa_freq_plot_output = gr.Image(label="Amino Acid Frequency Plot", type="pil", show_label=True, visible=False)
                seq_features_plot_output = gr.Image(label="Sequence Features Plot", type="pil", show_label=True, visible=False)
                seq_features_message_output = gr.Markdown() 
        with gr.TabItem("Pathways & Interactions"): # Renamed Tab
            with gr.Column(): # Using a single column for now for simplicity, can be changed to Row
                pathways_output = gr.Markdown(label="Biological Pathways (KEGG, Reactome)") # Pathways first
                interactions_output = gr.Markdown(label="Protein Interactions") # Then interactions

    submit_button.click(
        fn=get_protein_info,
        inputs=protein_id_input,
        outputs=[overview_output, aa_freq_plot_output, seq_features_plot_output, seq_features_message_output, pathways_output, interactions_output]
    )
    
    gr.Examples(
        examples=[["P05067"], ["P00533"], ["Q9BYF1"], ["P0DP23"], ["P04637"]],
        inputs=protein_id_input
    )

if __name__ == "__main__":
    iface.launch()
