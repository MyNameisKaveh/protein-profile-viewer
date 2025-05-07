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
        if feature['type'] not in plotted_feature_types_in_legend: # Use the original (case-sensitive) type for legend keys
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

# ------------------- START OF FINAL CORRECTED extract_sequence_features -------------------
def extract_sequence_features(uniprot_data, debug_markdown_list_ref): # debug_markdown_list_ref can be removed if not needed
    """
    Extracts and categorizes sequence features from UniProt data.
    Corrected for case-insensitivity and robust location parsing.
    """
    # features_of_interest keys are UPPERCASE for case-insensitive matching
    features_of_interest_uppercase = {
        "DOMAIN": "blue", "MOTIF": "green", "ACTIVE_SITE": "red",
        "BINDING_SITE": "orange", "MOD_RES": "purple", 
        "HELIX": "cyan", "STRAND": "magenta", "TURN": "gold"
    }
    extracted_features = []
    
    # Optional: Add a header for debug if you keep the debug list
    # debug_markdown_list_ref.append("\n\n--- FEATURE EXTRACTION (v5 - Case Corrected) ---\n")

    if "features" in uniprot_data and uniprot_data["features"] is not None:
        for feature_item in uniprot_data["features"]:
            feature_type_raw = feature_item.get("type")
            
            if not isinstance(feature_type_raw, str): # Skip if type is not a string
                continue

            feature_type_normalized = feature_type_raw.strip().upper() # Normalize to UPPERCASE and strip whitespace

            if feature_type_normalized in features_of_interest_uppercase:
                try:
                    location_obj = feature_item.get("location", {})
                    begin_pos_val_str = None; end_pos_val_str = None

                    start_node = location_obj.get("start")
                    end_node = location_obj.get("end")
                    position_node = location_obj.get("position")

                    if start_node and isinstance(start_node, dict) and "value" in start_node:
                        begin_pos_val_str = str(start_node["value"])
                    if end_node and isinstance(end_node, dict) and "value" in end_node:
                        end_pos_val_str = str(end_node["value"])
                    
                    if position_node and isinstance(position_node, dict) and "value" in position_node:
                        pos_val_str = str(position_node["value"])
                        if begin_pos_val_str is None: begin_pos_val_str = pos_val_str
                        if end_pos_val_str is None: end_pos_val_str = pos_val_str # Use for end if it's missing
                        # If it was truly a single position, make start and end the same
                        if "start" not in location_obj and "end" not in location_obj:
                             begin_pos_val_str = pos_val_str
                             end_pos_val_str = pos_val_str
                    
                    if begin_pos_val_str is None or end_pos_val_str is None:
                        continue

                    begin_pos = int(begin_pos_val_str)
                    end_pos = int(end_pos_val_str)

                    if begin_pos > end_pos:
                        continue
                        
                    extracted_features.append({
                        "type": feature_type_raw, # Store the original case for display/legend if needed
                        "begin": begin_pos, "end": end_pos,
                        "description": feature_item.get("description", feature_type_raw), 
                        "color": features_of_interest_uppercase[feature_type_normalized]
                    })
                except (ValueError, TypeError, AttributeError):
                    continue 
    return extracted_features
# ------------------- END OF FINAL CORRECTED extract_sequence_features -------------------

def get_protein_info(uniprot_id):
    if not uniprot_id:
        return "Please enter a UniProt ID.", gr.update(value=None, visible=False), gr.update(value=None, visible=False)
    url = UNIPROT_API_URL.format(accession=uniprot_id.strip().upper())
    try:
        response = requests.get(url); response.raise_for_status(); uniprot_api_data = response.json()
        primary_accession = uniprot_api_data.get("primaryAccession", "N/A")
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
        for comment in comments:
            if comment.get("commentType") == "FUNCTION":
                function_texts = comment.get("texts", [])
                if function_texts: function_comment_en = function_texts[0].get("value", "N/A"); break
        
        markdown_parts = [ # Using a list for easier construction
            f"**Primary Accession:** {primary_accession}\n", f"**Protein Name:** {protein_name}\n",
            f"**Gene Name(s):** {gene_name_str}\n", f"**Organism:** {organism_name}\n",
            f"**Sequence Length:** {length} amino acids\n", f"**Molecular Weight (calculated):** {calculated_mol_weight_str}\n\n",
            f"**Function Snippet:**\n{function_comment_en}\n\n",
            f"**Amino Acid Sequence (first 100 residues):**\n`{sequence[:100]}{'...' if len(sequence) > 100 else ''}`"
        ]
        
        # No more debug output to markdown_parts from here, extract_sequence_features will do its job
        # debug_markdown_list_ref can be removed from extract_sequence_features if no longer needed
        # For now, we pass an empty list if we don't want that specific debug output in the UI
        temp_debug_list = [] # Or remove this and the argument from extract_sequence_features
        seq_features = extract_sequence_features(uniprot_api_data, temp_debug_list) 
        
        # If you want to see the debug output in the UI, uncomment these lines:
        # if temp_debug_list:
        #     markdown_parts.extend(temp_debug_list)
        # markdown_parts.append(f"\n\n**Final Extracted features of interest count:** {len(seq_features)}\n")
        # if seq_features:
        #     extracted_feature_types_summary_counts = Counter([f['type'] for f in seq_features]) # Use original type for counting
        #     markdown_parts.append(f"**Final Extracted feature types (of interest) counts:** {dict(extracted_feature_types_summary_counts)}\n")

        aa_frequencies, freq_error = get_amino_acid_frequencies(sequence)
        aa_plot_update = gr.update(value=None, visible=False)
        if freq_error: markdown_parts.append(f"\n\n**Amino Acid Frequency Analysis:** {freq_error}")
        elif aa_frequencies:
            pil_image_aa = plot_amino_acid_frequencies(aa_frequencies)
            if pil_image_aa: aa_plot_update = gr.update(value=pil_image_aa, visible=True)
        
        feature_plot_update = gr.update(value=None, visible=False)
        if seq_features and length > 0:
            pil_image_features = plot_sequence_features(length, seq_features)
            if pil_image_features: feature_plot_update = gr.update(value=pil_image_features, visible=True)
            else: markdown_parts.append("\n\n**Sequence Features:** Could not generate feature plot.")
        elif not seq_features and length > 0 : markdown_parts.append("\n\n**Sequence Features:** No features of the selected types found or feature data is unavailable.")
        
        final_markdown = "".join(markdown_parts)
        return final_markdown, aa_plot_update, feature_plot_update

    except requests.exceptions.HTTPError as http_err: 
        status_code = http_err.response.status_code if http_err.response is not None else "N/A"
        error_message_en = f"Error: Protein with ID '{uniprot_id}' not found. Please check the ID." if status_code == 404 else f"HTTP error occurred: {http_err} (Status code: {status_code})"
        return error_message_en, gr.update(value=None, visible=False), gr.update(value=None, visible=False)
    except requests.exceptions.RequestException as req_err:
        return f"A network error occurred: {req_err}", gr.update(value=None, visible=False), gr.update(value=None, visible=False)
    except Exception as e:
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
    title="Protein Profile Viewer (v0.5 - Feature Fix Attempt)", 
    description="Enter a UniProt ID to display its basic information, amino acid frequency, and sequence features plot.",
    examples=[["P05067"], ["INS_HUMAN"], ["CYC_HUMAN"],["Q9BYF1"], ["P0DP23"], ["P00533"], ["P04637"]],
    flagging_options=None 
)
if __name__ == "__main__":
    iface.launch()
