import gradio as gr
import requests
from Bio.SeqUtils import molecular_weight
from collections import Counter
import matplotlib.pyplot as plt
import io
from PIL import Image
import sys 
import traceback
import json # For pretty printing JSON in debug

UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"

AMINO_ACID_NAMES = {
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic acid',
    'C': 'Cysteine', 'Q': 'Glutamine', 'E': 'Glutamic acid', 'G': 'Glycine',
    'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
    'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
    'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine'
}
STANDARD_AMINO_ACIDS_ORDER = "ARNDCQEGHILKMFPSTWYV" 

# --- Keep other functions (get_amino_acid_frequencies, plot_amino_acid_frequencies, plot_sequence_features) as they were in the last complete code ---
# ... (این توابع از کد کامل قبلی کپی شوند، بدون تغییر) ...

def get_amino_acid_frequencies(sequence):
    if not sequence or sequence == "N/A":
        return None, "Sequence not available for analysis."
    cleaned_sequence = "".join(filter(lambda x: x in AMINO_ACID_NAMES, sequence.upper()))
    if not cleaned_sequence:
        return None, "No valid amino acids found in sequence for counting."
    counts = Counter(cleaned_sequence)
    frequencies = {aa: counts.get(aa, 0) for aa in STANDARD_AMINO_ACIDS_ORDER}
    return frequencies, None

def plot_amino_acid_frequencies(frequencies):
    if not frequencies:
        return None
    ordered_keys = [key for key in STANDARD_AMINO_ACIDS_ORDER if key in frequencies]
    labels = [f"{aa}: {AMINO_ACID_NAMES.get(aa, aa)}" for aa in ordered_keys]
    values = [frequencies[aa] for aa in ordered_keys]
    fig, ax = plt.subplots(figsize=(12, 7)) 
    ax.bar(labels, values, color='skyblue')
    ax.set_xlabel("Amino Acid")
    ax.set_ylabel("Frequency")
    ax.set_title("Amino Acid Frequency Plot")
    plt.xticks(rotation=75, ha="right", fontsize=8) 
    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    img = Image.open(buf)
    plt.close(fig)
    return img

def plot_sequence_features(sequence_length, features):
    if not features or sequence_length == 0:
        return None
    fig, ax = plt.subplots(figsize=(12, max(3, len(features) * 0.4) + 1.5)) 
    ax.set_xlim(0, sequence_length)
    ax.set_xlabel("Amino Acid Position")
    ax.set_yticks([]) 
    ax.set_title("Sequence Features Plot")
    legend_handles = {} 
    y_pos_counter = 0
    plotted_feature_types_in_legend = set()
    bar_height = 0.8 
    for feature in sorted(features, key=lambda x: x["begin"]):
        begin = feature["begin"]
        end = feature["end"]
        color = feature["color"]
        width = max(1, end - begin + 1) 
        ax.barh(y_pos_counter, width, height=bar_height, left=begin -1, color=color, edgecolor='black', alpha=0.7)
        if feature['type'] not in plotted_feature_types_in_legend:
            legend_handles[feature['type']] = plt.Rectangle((0, 0), 1, 1, fc=color, alpha=0.7)
            plotted_feature_types_in_legend.add(feature['type'])
        y_pos_counter += 1 
    if y_pos_counter > 0:
        ax.set_ylim(-0.5, y_pos_counter -1 + bar_height/2 + 0.5) 
    else: 
        plt.close(fig)
        return None
    if legend_handles:
        ax.legend(legend_handles.values(), legend_handles.keys(), title="Feature Types", bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    plt.tight_layout(rect=[0, 0, 0.83, 0.96])
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    img = Image.open(buf)
    plt.close(fig)
    return img
# --- END OF UNCHANGED FUNCTIONS ---


# ------------------- START OF HIGHLY DEBUGGABLE extract_sequence_features -------------------
def extract_sequence_features(uniprot_data, debug_markdown_list_ref):
    """
    Extracts and categorizes sequence features from UniProt data.
    Includes extensive debugging output.
    """
    features_of_interest = {
        "DOMAIN": "blue", "MOTIF": "green", "ACTIVE_SITE": "red",
        "BINDING_SITE": "orange", "MOD_RES": "purple", 
        "HELIX": "cyan", "STRAND": "magenta", "TURN": "gold"
    }
    
    extracted_features = []
    # Add debug header to the list that will be joined into markdown
    debug_markdown_list_ref.append("\n\n--- DETAILED FEATURE EXTRACTION DEBUG ---\n")

    if "features" in uniprot_data and uniprot_data["features"] is not None:
        # Limit debugging to a few features to avoid excessively long output
        features_to_debug_count = 0
        max_features_to_log = 5 # Log details for first 5 relevant features

        for i, feature in enumerate(uniprot_data["features"]):
            feature_type = feature.get("type")
            location_obj = feature.get("location", {}) # Get location object early

            # Log only if it's a type we are interested in, or for the first few overall
            if feature_type in features_of_interest and features_to_debug_count < max_features_to_log:
                features_to_debug_count += 1
                debug_markdown_list_ref.append(f"\n**Processing Raw Feature Index {i} (Type: {feature_type}):**\n")
                debug_markdown_list_ref.append(f"```json\nLocation Object: {json.dumps(location_obj, indent=2)}\n```\n")
                
                begin_pos_val_str = None
                end_pos_val_str = None

                # Attempt 1: Standard 'start' and 'end' with 'value'
                start_node = location_obj.get("start")
                end_node = location_obj.get("end")
                
                if start_node is not None and isinstance(start_node, dict) and "value" in start_node:
                    begin_pos_val_str = str(start_node["value"])
                    debug_markdown_list_ref.append(f"- Attempt 1 (start.value): `begin_pos_val_str` = `{begin_pos_val_str}`\n")
                else:
                    debug_markdown_list_ref.append(f"- Attempt 1 (start.value): Not found or invalid structure. `start_node`: {start_node}\n")

                if end_node is not None and isinstance(end_node, dict) and "value" in end_node:
                    end_pos_val_str = str(end_node["value"])
                    debug_markdown_list_ref.append(f"- Attempt 1 (end.value): `end_pos_val_str` = `{end_pos_val_str}`\n")
                else:
                    debug_markdown_list_ref.append(f"- Attempt 1 (end.value): Not found or invalid structure. `end_node`: {end_node}\n")

                # Attempt 2: Single 'position' with 'value'
                if begin_pos_val_str is None and end_pos_val_str is None and \
                   "position" in location_obj and isinstance(location_obj["position"], dict) and \
                   "value" in location_obj["position"]:
                    pos_val_str = str(location_obj["position"]["value"])
                    begin_pos_val_str = pos_val_str
                    end_pos_val_str = pos_val_str
                    debug_markdown_list_ref.append(f"- Attempt 2 (position.value): `begin_pos_val_str` = `{begin_pos_val_str}`, `end_pos_val_str` = `{end_pos_val_str}`\n")
                elif begin_pos_val_str is None and end_pos_val_str is None: # Only if attempt 1 failed completely
                     debug_markdown_list_ref.append(f"- Attempt 2 (position.value): Not found or invalid structure. `location_obj.get('position')`: {location_obj.get('position')}\n")


                if begin_pos_val_str is None or end_pos_val_str is None:
                    debug_markdown_list_ref.append(f"-> **Result for Feature Index {i}: SKIPPED (Could not determine valid begin/end)**\n")
                    continue

                try:
                    begin_pos = int(begin_pos_val_str)
                    end_pos = int(end_pos_val_str)
                    debug_markdown_list_ref.append(f"- Parsed Integers: `begin_pos` = {begin_pos}, `end_pos` = {end_pos}\n")

                    if begin_pos > end_pos:
                        debug_markdown_list_ref.append(f"-> **Result for Feature Index {i}: SKIPPED (begin > end: {begin_pos} > {end_pos})**\n")
                        continue

                    description = feature.get("description", feature_type)
                    extracted_features.append({
                        "type": feature_type, "begin": begin_pos, "end": end_pos,
                        "description": description, "color": features_of_interest[feature_type]
                    })
                    debug_markdown_list_ref.append(f"-> **Result for Feature Index {i}: ADDED TO EXTRACTED LIST**\n")

                except (ValueError, TypeError) as e_conv: # Catch conversion or type errors
                    debug_markdown_list_ref.append(f"-> **Result for Feature Index {i}: SKIPPED (Error during int conversion: {e_conv})**\n")
                    continue
            
            elif feature_type in features_of_interest: # It was of interest but we already logged max_features_to_log
                # Still try to extract it, just don't log verbosely
                try:
                    location_obj_silent = feature.get("location", {})
                    begin_pos_val_str_silent = None; end_pos_val_str_silent = None
                    start_node_silent = location_obj_silent.get("start"); end_node_silent = location_obj_silent.get("end")
                    if start_node_silent and "value" in start_node_silent: begin_pos_val_str_silent = str(start_node_silent["value"])
                    if end_node_silent and "value" in end_node_silent: end_pos_val_str_silent = str(end_node_silent["value"])
                    if begin_pos_val_str_silent is None and end_pos_val_str_silent is None and "position" in location_obj_silent and "value" in location_obj_silent["position"]:
                        pos_val_str_silent = str(location_obj_silent["position"]["value"])
                        begin_pos_val_str_silent = pos_val_str_silent; end_pos_val_str_silent = pos_val_str_silent
                    if begin_pos_val_str_silent and end_pos_val_str_silent:
                        bp_s = int(begin_pos_val_str_silent); ep_s = int(end_pos_val_str_silent)
                        if bp_s <= ep_s:
                            extracted_features.append({
                                "type": feature_type, "begin": bp_s, "end": ep_s,
                                "description": feature.get("description", feature_type), 
                                "color": features_of_interest[feature_type]
                            })
                except: continue # Silent fail for non-logged features

    debug_markdown_list_ref.append("\n--- END OF DETAILED FEATURE EXTRACTION DEBUG ---\n")
    return extracted_features
# ------------------- END OF HIGHLY DEBUGGABLE extract_sequence_features -------------------


def get_protein_info(uniprot_id):
    if not uniprot_id:
        return "Please enter a UniProt ID.", gr.update(value=None, visible=False), gr.update(value=None, visible=False)
    url = UNIPROT_API_URL.format(accession=uniprot_id.strip().upper())
    try:
        response = requests.get(url)
        response.raise_for_status()
        uniprot_api_data = response.json()

        primary_accession = uniprot_api_data.get("primaryAccession", "N/A")
        protein_data_dict = uniprot_api_data.get("proteinDescription", {}).get("recommendedName", {})
        protein_name = protein_data_dict.get("fullName", {}).get("value", "N/A")
        if protein_name == "N/A" and uniprot_api_data.get("proteinDescription", {}).get("submissionNames"):
            protein_name = uniprot_api_data["proteinDescription"]["submissionNames"][0].get("fullName", {}).get("value", "N/A")
        gene_names_data = uniprot_api_data.get("genes")
        gene_name_str = "N/A"
        if gene_names_data:
            gene_names_list = [g.get("geneName", {}).get("value", "") for g in gene_names_data if g.get("geneName")]
            if not gene_names_list:
                 gene_names_list = [g.get("orfNames", [{}])[0].get("value", "") for g in gene_names_data if g.get("orfNames")]
            gene_name_str = ", ".join(filter(None, gene_names_list)) or "N/A"
        organism_data = uniprot_api_data.get("organism", {})
        organism_name = organism_data.get("scientificName", "N/A")
        if organism_data.get("commonName"):
            organism_name += f" ({organism_data.get('commonName')})"
        sequence_info = uniprot_api_data.get("sequence", {})
        sequence = sequence_info.get("value", "N/A")
        length = sequence_info.get("length", 0)
        calculated_mol_weight_str = "N/A"
        if sequence != "N/A" and length > 0:
            try:
                cleaned_sequence_for_mw = "".join(filter(lambda x: x in AMINO_ACID_NAMES, sequence.upper()))
                if cleaned_sequence_for_mw:
                     calculated_mol_weight = molecular_weight(cleaned_sequence_for_mw, seq_type="protein")
                     calculated_mol_weight_str = f"{calculated_mol_weight:.2f} Da"
                else:
                    calculated_mol_weight_str = "Invalid sequence for MW calculation"
            except Exception:
                calculated_mol_weight_str = "Error in MW calculation"
        comments = uniprot_api_data.get("comments", [])
        function_comment_en = "N/A" 
        for comment in comments:
            if comment.get("commentType") == "FUNCTION":
                function_texts = comment.get("texts", [])
                if function_texts:
                    function_comment_en = function_texts[0].get("value", "N/A")
                    break
        
        # Initialize a list to gather all markdown parts, including debug info
        markdown_parts = [
            f"**Primary Accession:** {primary_accession}\n",
            f"**Protein Name:** {protein_name}\n",
            f"**Gene Name(s):** {gene_name_str}\n",
            f"**Organism:** {organism_name}\n",
            f"**Sequence Length:** {length} amino acids\n",
            f"**Molecular Weight (calculated):** {calculated_mol_weight_str}\n\n",
            f"**Function Snippet:**\n{function_comment_en}\n\n",
            f"**Amino Acid Sequence (first 100 residues):**\n`{sequence[:100]}{'...' if len(sequence) > 100 else ''}`"
        ]
        
        # ---- GENERAL DEBUG INFO ----
        markdown_parts.append(f"\n\n--- GENERAL DEBUG INFO ---\n")
        markdown_parts.append(f"**Raw features from API count:** {len(uniprot_api_data.get('features', []))}\n")
        if uniprot_api_data.get('features'):
            raw_feature_types = Counter([f.get('type', 'UNKNOWN_TYPE') for f in uniprot_api_data['features']])
            markdown_parts.append(f"**Raw feature types found:** {dict(raw_feature_types)}\n")
        
        # Pass the list to extract_sequence_features so it can append detailed debug info
        seq_features = extract_sequence_features(uniprot_api_data, markdown_parts) 
        
        markdown_parts.append(f"**Extracted features (of interest) count (after detailed debug):** {len(seq_features)}\n")
        if seq_features:
            extracted_feature_types_counts = Counter([f['type'] for f in seq_features])
            markdown_parts.append(f"**Extracted feature types (of interest) counts:** {dict(extracted_feature_types_counts)}\n")
        markdown_parts.append(f"--- END GENERAL DEBUG INFO ---\n")
        # ---- END GENERAL DEBUG INFO ----
        
        # Amino acid frequency plot
        aa_frequencies, freq_error = get_amino_acid_frequencies(sequence)
        aa_plot_update = gr.update(value=None, visible=False)
        if freq_error:
            markdown_parts.append(f"\n\n**Amino Acid Frequency Analysis:** {freq_error}")
        elif aa_frequencies:
            pil_image_aa = plot_amino_acid_frequencies(aa_frequencies)
            if pil_image_aa:
                aa_plot_update = gr.update(value=pil_image_aa, visible=True)

        # Sequence features plot
        feature_plot_update = gr.update(value=None, visible=False)
        if seq_features and length > 0:
            pil_image_features = plot_sequence_features(length, seq_features)
            if pil_image_features:
                feature_plot_update = gr.update(value=pil_image_features, visible=True)
            else:
                 markdown_parts.append("\n\n**Sequence Features:** Could not generate feature plot (e.g. no plottable features after filtering).")
        elif not seq_features and length > 0 : 
             markdown_parts.append("\n\n**Sequence Features:** No features of the selected types found or feature data is unavailable.")
        
        final_markdown = "".join(markdown_parts)
        return final_markdown, aa_plot_update, feature_plot_update

    except requests.exceptions.HTTPError as http_err:
        status_code = http_err.response.status_code if http_err.response is not None else "N/A"
        error_message_en = f"Error: Protein with ID '{uniprot_id}' not found. Please check the ID." if status_code == 404 else f"HTTP error occurred: {http_err} (Status code: {status_code})"
        return error_message_en, gr.update(value=None, visible=False), gr.update(value=None, visible=False)
    except requests.exceptions.RequestException as req_err:
        return f"A network error occurred: {req_err}", gr.update(value=None, visible=False), gr.update(value=None, visible=False)
    except Exception as e:
        # print(f"Unexpected error for {uniprot_id}: {e}", file=sys.stderr) 
        # traceback.print_exc(file=sys.stderr) 
        return f"An unexpected error occurred while processing ID '{uniprot_id}'. Details: {str(e)[:150]}", gr.update(value=None, visible=False), gr.update(value=None, visible=False)

outputs_list = [
    gr.Markdown(label="Protein Information"),
    gr.Image(label="Amino Acid Frequency Plot", type="pil", show_label=True, visible=False),
    gr.Image(label="Sequence Features Plot", type="pil", show_label=True, visible=False)
]

iface = gr.Interface(
    fn=get_protein_info,
    inputs=gr.Textbox(label="Enter UniProt ID (e.g., P05067 or INS_HUMAN)", placeholder="e.g., P0DP23"),
    outputs=outputs_list,
    title="Protein Profile Viewer (v0.4.3 - Extensive Debug)", # Version updated
    description="Enter a UniProt ID to display its basic information, amino acid frequency, and sequence features plot. Extensive debug info for features is currently active.",
    examples=[["P05067"], ["INS_HUMAN"], ["CYC_HUMAN"],["Q9BYF1"], ["P0DP23"], ["P00533"], ["P04637"]],
    flagging_options=None 
)

if __name__ == "__main__":
    iface.launch()
