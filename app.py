import gradio as gr
import requests
from Bio.SeqUtils import molecular_weight
from collections import Counter
import matplotlib.pyplot as plt
import io
from PIL import Image
import sys # For more detailed error logging if needed
import traceback # For more detailed error logging if needed

UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"

def get_amino_acid_frequencies(sequence):
    """
    Calculates the frequency of each amino acid in the sequence.
    """
    if not sequence or sequence == "N/A":
        return None, "Sequence not available for analysis."
    
    standard_amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    cleaned_sequence = "".join(filter(lambda x: x in standard_amino_acids, sequence.upper()))
    
    if not cleaned_sequence:
        return None, "No valid amino acids found in sequence for counting."
        
    counts = Counter(cleaned_sequence)
    frequencies = {aa: counts.get(aa, 0) for aa in standard_amino_acids}
    return frequencies, None

def plot_amino_acid_frequencies(frequencies):
    """
    Plots a bar chart of amino acid frequencies and returns a PIL Image object.
    """
    if not frequencies:
        return None

    names = list(frequencies.keys())
    values = list(frequencies.values())

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(names, values, color='skyblue')
    ax.set_xlabel("Amino Acid")
    ax.set_ylabel("Frequency")
    ax.set_title("Amino Acid Frequency Plot")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    img = Image.open(buf)
    plt.close(fig)
    return img


def get_protein_info(uniprot_id):
    """
    Fetches and processes protein information from UniProt API.
    Returns main info markdown and a PIL image for the frequency plot (or updates for Gradio components).
    """
    if not uniprot_id:
        # Return appropriate number of Nones or empty updates for all outputs
        return "Please enter a UniProt ID.", gr.update(value=None, visible=False)

    url = UNIPROT_API_URL.format(accession=uniprot_id.strip().upper())

    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()

        primary_accession = data.get("primaryAccession", "N/A")
        protein_data_dict = data.get("proteinDescription", {}).get("recommendedName", {})
        protein_name = protein_data_dict.get("fullName", {}).get("value", "N/A")
        if protein_name == "N/A" and data.get("proteinDescription", {}).get("submissionNames"):
            protein_name = data["proteinDescription"]["submissionNames"][0].get("fullName", {}).get("value", "N/A")

        gene_names_data = data.get("genes")
        gene_name_str = "N/A"
        if gene_names_data:
            gene_names_list = [g.get("geneName", {}).get("value", "") for g in gene_names_data if g.get("geneName")]
            if not gene_names_list:
                 gene_names_list = [g.get("orfNames", [{}])[0].get("value", "") for g in gene_names_data if g.get("orfNames")]
            gene_name_str = ", ".join(filter(None, gene_names_list)) or "N/A"


        organism_data = data.get("organism", {})
        organism_name = organism_data.get("scientificName", "N/A")
        if organism_data.get("commonName"):
            organism_name += f" ({organism_data.get('commonName')})"

        sequence_info = data.get("sequence", {})
        sequence = sequence_info.get("value", "N/A")
        length = sequence_info.get("length", 0)
        
        calculated_mol_weight_str = "N/A"
        if sequence != "N/A":
            try:
                cleaned_sequence_for_mw = "".join(filter(lambda x: x in "ACDEFGHIKLMNPQRSTVWY", sequence.upper()))
                if cleaned_sequence_for_mw:
                     calculated_mol_weight = molecular_weight(cleaned_sequence_for_mw, seq_type="protein")
                     calculated_mol_weight_str = f"{calculated_mol_weight:.2f} Da"
                else:
                    calculated_mol_weight_str = "Invalid sequence for MW calculation"
            except Exception:
                calculated_mol_weight_str = "Error in MW calculation"

        comments = data.get("comments", [])
        function_comment_en = "N/A" # Changed to English
        for comment in comments:
            if comment.get("commentType") == "FUNCTION":
                function_texts = comment.get("texts", [])
                if function_texts:
                    function_comment_en = function_texts[0].get("value", "N/A")
                    break
        
        # Main protein information in Markdown (English)
        main_info_md = (
            f"**Primary Accession:** {primary_accession}\n"
            f"**Protein Name:** {protein_name}\n"
            f"**Gene Name(s):** {gene_name_str}\n"
            f"**Organism:** {organism_name}\n"
            f"**Sequence Length:** {length} amino acids\n"
            f"**Molecular Weight (calculated):** {calculated_mol_weight_str}\n\n"
            f"**Function Snippet:**\n{function_comment_en}\n\n"
            f"**Amino Acid Sequence (first 100 residues):**\n`{sequence[:100]}{'...' if len(sequence) > 100 else ''}`"
        )
        
        aa_frequencies, freq_error = get_amino_acid_frequencies(sequence)
        aa_plot_pil_image_update = gr.update(value=None, visible=False) # Default to hidden

        if freq_error:
            main_info_md += f"\n\n**Amino Acid Frequency Analysis:** {freq_error}"
        elif aa_frequencies:
            pil_image = plot_amino_acid_frequencies(aa_frequencies)
            if pil_image:
                aa_plot_pil_image_update = gr.update(value=pil_image, visible=True)
        
        return main_info_md, aa_plot_pil_image_update

    except requests.exceptions.HTTPError as http_err:
        status_code = http_err.response.status_code if http_err.response is not None else "N/A"
        error_message_en = f"Error: Protein with ID '{uniprot_id}' not found. Please check the ID." if status_code == 404 else f"HTTP error occurred: {http_err} (Status code: {status_code})"
        return error_message_en, gr.update(value=None, visible=False)
    except requests.exceptions.RequestException as req_err:
        return f"A network error occurred: {req_err}", gr.update(value=None, visible=False)
    except Exception as e:
        # For debugging, you might want to log the full error
        # print(f"Unexpected error for {uniprot_id}: {e}", file=sys.stderr)
        # traceback.print_exc(file=sys.stderr)
        return f"An unexpected error occurred while processing the data. Please check the ID and try again.", gr.update(value=None, visible=False)


# Define Gradio UI (English)
outputs_list = [
    gr.Markdown(label="Protein Information"),
    gr.Image(label="Amino Acid Frequency Plot", type="pil", show_label=True, visible=False) # Start hidden
]

iface = gr.Interface(
    fn=get_protein_info,
    inputs=gr.Textbox(label="Enter UniProt ID (e.g., P05067 or INS_HUMAN)", placeholder="e.g., P0DP23"),
    outputs=outputs_list,
    title="Protein Profile Viewer (Enhanced)",
    description="Enter a UniProt ID to display its basic information and amino acid frequency plot.",
    examples=[["P05067"], ["INS_HUMAN"], ["CYC_HUMAN"], ["P12345"]],
    flagging_options=None 
)

if __name__ == "__main__":
    iface.launch()
