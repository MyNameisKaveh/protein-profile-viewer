# -*- coding: utf-8 -*-
# Import necessary libraries
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
import pandas as pd # Keep import, might be useful later

# --- Constants ---
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
    """Calculates the frequency of each standard amino acid in the sequence."""
    if not sequence or sequence == "N/A": return None, "Sequence not available for analysis."
    # Filter sequence for standard amino acids and count
    cleaned_sequence = "".join(filter(lambda x: x in AMINO_ACID_NAMES, sequence.upper()))
    if not cleaned_sequence: return None, "No valid standard amino acids found for counting."
    counts = Counter(cleaned_sequence)
    # Ensure all 20 standard AAs are in the dict with defined order
    frequencies = {aa: counts.get(aa, 0) for aa in STANDARD_AMINO_ACIDS_ORDER}
    return frequencies, None

def plot_amino_acid_frequencies(frequencies):
    """Plots a bar chart of amino acid frequencies."""
    if not frequencies: return None
    ordered_keys = [key for key in STANDARD_AMINO_ACIDS_ORDER if key in frequencies]
    labels = [f"{aa}: {AMINO_ACID_NAMES.get(aa, aa)}" for aa in ordered_keys]
    values = [frequencies[aa] for aa in ordered_keys]
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 7)) 
    ax.bar(labels, values, color='skyblue')
    ax.set_xlabel("Amino Acid")
    ax.set_ylabel("Frequency")
    ax.set_title("Amino Acid Frequency Plot")
    plt.xticks(rotation=75, ha="right", fontsize=8) 
    plt.tight_layout() # Adjust layout
    
    # Save plot to buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    img = Image.open(buf)
    plt.close(fig) # Close plot to free memory
    return img

def extract_sequence_features(uniprot_data):
    """Extracts and categorizes sequence features."""
    features_of_interest_uppercase = { # Use uppercase keys for matching
        "DOMAIN": "blue", "MOTIF": "green", "ACTIVE_SITE": "red",
        "BINDING_SITE": "orange", "MOD_RES": "purple", 
        "HELIX": "cyan", "STRAND": "magenta", "TURN": "gold"
    }
    extracted = []
    if "features" in uniprot_data and uniprot_data["features"]:
        for item in uniprot_data["features"]:
            type_raw = item.get("type")
            if not isinstance(type_raw, str): continue
            type_norm = type_raw.strip().upper() # Normalize for matching
            
            if type_norm in features_of_interest_uppercase:
                try:
                    loc_obj = item.get("location", {})
                    b_str, e_str = None, None
                    s_node, e_node, p_node = loc_obj.get("start"), loc_obj.get("end"), loc_obj.get("position")
                    
                    # Robust location parsing
                    if s_node and isinstance(s_node, dict) and "value" in s_node: b_str = str(s_node["value"])
                    if e_node and isinstance(e_node, dict) and "value" in e_node: e_str = str(e_node["value"])
                    if p_node and isinstance(p_node, dict) and "value" in p_node:
                        p_str = str(p_node["value"])
                        if b_str is None: b_str = p_str
                        if e_str is None: e_str = p_str
                        if "start" not in loc_obj and "end" not in loc_obj: b_str, e_str = p_str, p_str
                            
                    if b_str is None or e_str is None: continue # Skip if cannot determine position
                    
                    b_pos, e_pos = int(b_str), int(e_str)
                    if b_pos > e_pos: continue # Skip invalid range
                    
                    extracted.append({
                        "type": type_raw, # Store original type for display
                        "begin": b_pos, "end": e_pos, 
                        "description": item.get("description", type_raw), 
                        "color": features_of_interest_uppercase[type_norm] # Get color using normalized type
                    })
                except (ValueError, TypeError, AttributeError): 
                    continue # Skip on parsing errors
    return extracted

def plot_sequence_features(sequence_length, features):
    """Plots sequence features on a horizontal bar."""
    if not features or sequence_length == 0: return None
    
    # Adjust figure size based on number of features
    fig, ax = plt.subplots(figsize=(12, max(3, len(features) * 0.4) + 1.5)) 
    ax.set_xlim(0, sequence_length)
    ax.set_xlabel("Amino Acid Position")
    ax.set_yticks([]) # Hide y-axis ticks
    ax.set_title("Sequence Features Plot")
    
    leg_h = {}; y_pos, leg_set, b_h = 0, set(), 0.8 
    
    # Plot features sorted by start position
    for feat in sorted(features, key=lambda x: x["begin"]):
        b, e, c = feat["begin"], feat["end"], feat["color"]
        w = max(1, e - b + 1) # Ensure minimum width of 1
        # Plot bar (left position is 0-based)
        ax.barh(y_pos, w, height=b_h, left=b - 1, color=c, edgecolor='black', alpha=0.7)
        
        # Add feature type to legend if new
        if feat['type'] not in leg_set: 
            leg_h[feat['type']] = plt.Rectangle((0, 0), 1, 1, fc=c, alpha=0.7)
            leg_set.add(feat['type'])
        y_pos += 1 # Move to next vertical position
        
    if y_pos > 0: 
        ax.set_ylim(-0.5, y_pos -1 + b_h/2 + 0.5) # Adjust y-limits
    else: 
        plt.close(fig); return None # No features plotted
        
    # Add legend outside the plot area
    if leg_h: 
        ax.legend(leg_h.values(), leg_h.keys(), title="Feature Types", 
                  bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
                  
    plt.tight_layout(rect=[0, 0, 0.83, 0.96]) # Adjust layout for legend
    
    # Save plot to buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    img = Image.open(buf)
    plt.close(fig) # Close plot to free memory
    return img

def extract_interactions(uniprot_data):
    """Extracts protein interaction partners from comments."""
    interactions = []
    if "comments" in uniprot_data:
        for c in uniprot_data["comments"]:
            if c.get("commentType") == "INTERACTION" and "interactions" in c:
                for i_entry in c["interactions"]:
                    i1_acc = i_entry.get("interactantOne", {}).get("uniProtKBAccession")
                    i2_acc = i_entry.get("interactantTwo", {}).get("uniProtKBAccession")
                    i2_gene = i_entry.get("interactantTwo", {}).get("geneName")
                    # Check if the current protein is interactantOne and partner exists
                    if i1_acc == uniprot_data.get("primaryAccession") and i2_acc:
                        partner = f"{i2_gene} ({i2_acc})" if i2_gene else i2_acc
                        interactions.append(f"- Interacts with: **{partner}**")
    return "\n".join(sorted(interactions)) if interactions else "No interaction partners listed in UniProt comments."

def extract_pathways(uniprot_data):
    """Extracts pathway links from cross-references."""
    pathways = []
    # URLs for KEGG and Reactome (simplified for direct use)
    dbs = { 
        "KEGG": "https://www.genome.jp/dbget-bin/www_bget?", 
        "Reactome": "https://reactome.org/content/detail/" 
    }
    if "uniProtKBCrossReferences" in uniprot_data:
        for xref in uniprot_data["uniProtKBCrossReferences"]:
            db = xref.get("database")
            pid = xref.get("id")
            if db in dbs and pid:
                desc = pid # Default display text is the ID
                # Try to get a more descriptive name from properties
                if "properties" in xref:
                    for p in xref["properties"]:
                        if p.get("key") in ["PathwayName", "Description"]: 
                            desc = f"{p.get('value')} ({pid})"
                            break
                link = dbs[db] + pid # Construct the link
                pathways.append(f"- [{desc}]({link}) ({db})") # Format as markdown link
    return "\n".join(sorted(list(set(pathways)))) if pathways else "No KEGG/Reactome pathway info found."

def extract_disease_info(uniprot_data):
    """Extracts disease association information."""
    diseases = []
    if "comments" in uniprot_data:
        for c in uniprot_data["comments"]:
            if c.get("commentType") == "DISEASE" and "disease" in c:
                d_entry = c["disease"]; d_name = d_entry.get("diseaseId", "Unknown"); desc = d_entry.get("description", "N/A")
                mim = None
                # Extract MIM link if available
                if "diseaseCrossReference" in d_entry and d_entry["diseaseCrossReference"].get("database") == "MIM": 
                    mim = d_entry["diseaseCrossReference"].get("id")
                d_md = f"**{d_name}**" + (f" (MIM: [{mim}](https://omim.org/entry/{mim}))" if mim else "")
                d_md += f"\n   - *Desc:* {desc}\n"
                # Extract note if available
                note_val = c.get("note", {}).get("texts", [{}])[0].get("value")
                if note_val: d_md += f"   - *Note:* {note_val}\n"
                diseases.append(d_md)
    return "\n---\n".join(sorted(diseases)) if diseases else "No disease association info found."

def extract_publications(uniprot_data):
    """Extracts publication information."""
    pubs = []
    if "references" in uniprot_data:
        for i, ref in enumerate(uniprot_data.get("references", [])):
            cit = ref.get("citation", {}); title = cit.get("title", "N/A"); authors = ", ".join(cit.get("authors", ["N/A"]))
            j = cit.get("journalName", ""); v = cit.get("volume", ""); f = cit.get("firstPage", ""); l = cit.get("lastPage", ""); d = cit.get("publicationDate", "")
            pmid, doi = None, None
            if "citationCrossReferences" in cit:
                for xr in cit["citationCrossReferences"]:
                    if xr.get("database") == "PubMed": pmid = xr.get("id")
                    elif xr.get("database") == "DOI": doi = xr.get("id")
            
            # Format publication entry
            md = f"**{i + 1}. {title}**\n   - *{authors}*\n   - *{j}" + (f", {v}" if v else "") + (f":{f}" if f else "") + (f"-{l}" if l else "")
            # Correctly add date and newline
            if d: 
                md += f" ({d})" 
            md += "*\n" # End italics and add newline
            
            if pmid: md += f"   - [PubMed {pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)\n"
            if doi: md += f"   - [DOI {doi}](https://doi.org/{doi})\n"
            pubs.append(md)
    return "\n---\n".join(pubs) if pubs else "No publication info found."

def extract_cross_references(uniprot_data):
    """Extracts links to selected external databases."""
    xrefs = []; 
    # Define target databases and URL patterns (simplified)
    dbs = {
        "Ensembl": "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", 
        "GeneID": "https://www.ncbi.nlm.nih.gov/gene/", 
        "RefSeq": "https://www.ncbi.nlm.nih.gov/nuccore/", 
        "GO": "https://amigo.geneontology.org/amigo/term/", 
        "InterPro": "https://www.ebi.ac.uk/interpro/entry/InterPro/", 
        "Pfam": "https://www.ebi.ac.uk/interpro/entry/pfam/", 
        "PDB": "https://www.rcsb.org/structure/",
        "KEGG": "https://www.genome.jp/dbget-bin/www_bget?", 
        "Reactome": "https://reactome.org/content/detail/"
    }
    grouped = {db: [] for db in dbs}; # Group links by database
    
    if "uniProtKBCrossReferences" in uniprot_data:
        for xr in uniprot_data["uniProtKBCrossReferences"]:
            db = xr.get("database"); xid = xr.get("id")
            if db in dbs and xid:
                url, txt = None, xid # Default display text is ID
                # Specific logic for some DBs to get better links/text
                if db == "Ensembl": 
                    gene_id = next((p.get("value") for p in xr.get("properties", []) if p.get("key") == "GeneId"), None)
                    if gene_id: url, txt = dbs[db] + gene_id, gene_id
                elif db == "RefSeq":
                    prot_id = next((p.get("value") for p in xr.get("properties", []) if p.get("key") == "ProteinId"), None)
                    if prot_id: url, txt = dbs[db] + prot_id, prot_id
                elif db == "GO":
                    term = next((p.get("value") for p in xr.get("properties", []) if p.get("key") == "GoTerm"), xid)
                    url, txt = dbs[db] + xid, f"{term} ({xid})"
                else: # Default URL construction
                    url = dbs[db] + xid
                    
                if url: 
                    # Avoid duplicates within the same database category
                    link_md = f"[{txt}]({url})"
                    if link_md not in grouped[db]:
                         grouped[db].append(link_md)

        # Format the grouped links
        for db, links in grouped.items():
            if links: xrefs.append(f"**{db}:** " + ", ".join(sorted(links)))
            
    return "\n".join(xrefs) if xrefs else "No cross-references found for the selected databases."

# --- Search Function ---
def search_uniprot_by_name(search_term, result_limit=5):
    """Searches UniProtKB for a protein/gene name."""
    if not search_term or len(search_term) < 3: return "Enter at least 3 characters.", []
    params = { "query": f'({search_term}) AND (reviewed:true)', "fields": "accession,id,protein_name,organism_name", "format": "json", "size": result_limit }
    status_md = f"### Search Results for '{search_term}':\n"; dataset_data = []
    try:
        response = requests.get(UNIPROT_SEARCH_URL, params=params); response.raise_for_status(); data = response.json()
        results = data.get("results")
        if not results: status_md += "No reviewed entries found."
        else:
            status_md += f"*Found {len(results)}. Select row, then copy ID.*\n---"
            for entry in results:
                acc = entry.get("primaryAccession", "N/A"); uid = entry.get("uniProtkbId", ""); name = "N/A"
                if entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"): name = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
                elif entry.get("proteinDescription", {}).get("submissionNames"): name = entry["proteinDescription"]["submissionNames"][0].get("fullName",{}).get("value", "N/A")
                org = entry.get("organism", {}).get("scientificName", "N/A"); dataset_data.append([acc, uid, name, org])
            status_md += "\n---" 
    except Exception as e: status_md += f"\n**Search Error:** {e}"
    return status_md, dataset_data

# --- Function to update copy box on Dataset selection ---
def update_copy_box(evt: gr.SelectData):
    """Gets the Accession ID from the selected row."""
    # The SelectData event for Dataset has 'index' (row, col) and 'value' (cell content)
    # We assume the Accession ID is in the first column (index 0)
    # evt.index[1] should be 0 for the first column.
    # evt.value will contain the value of the selected cell (which is the Accession ID if first cell clicked)
    # A safer way might be to use evt.index[0] (row index) and get the data from the dataset component
    # but using evt.value is simpler if the user clicks the first column.
    if evt.value and evt.index[1] == 0: # Check if a value exists and it's from the first column
        return evt.value
    return "" # Return empty string otherwise

# --- Main Function to Get All Info ---
def get_protein_info(uniprot_id):
    """Fetches and processes all protein information for display."""
    empty_plot = gr.update(value=None, visible=False); empty_str = ""; err_msg = "Error"
    # Define the number of outputs expected
    num_outputs = 9 
    outputs_on_error = tuple([err_msg] + [empty_plot if i in [1, 2] else empty_str for i in range(1, num_outputs)])

    if not uniprot_id: 
        return ("Enter UniProt ID.",) + outputs_on_error[1:]
    
    url = UNIPROT_API_URL.format(accession=uniprot_id.strip().upper())
    try:
        response = requests.get(url); response.raise_for_status(); data = response.json()
        
        acc = data.get("primaryAccession", "N/A"); id_display = data.get("uniProtkbId", "N/A")
        link = f"https://www.uniprot.org/uniprotkb/{acc}/entry"; name_dict = data.get("proteinDescription", {}).get("recommendedName", {}); name = name_dict.get("fullName", {}).get("value", "N/A")
        if name == "N/A" and data.get("proteinDescription", {}).get("submissionNames"): name = data["proteinDescription"]["submissionNames"][0].get("fullName", {}).get("value", "N/A")
        genes_data = data.get("genes"); gene_str = "N/A"
        if genes_data:
            g_list = [g.get("geneName", {}).get("value", "") for g in genes_data if g.get("geneName")]
            if not g_list: g_list = [g.get("orfNames", [{}])[0].get("value", "") for g in genes_data if g.get("orfNames")]
            gene_str = ", ".join(filter(None, g_list)) or "N/A"
        org_data = data.get("organism", {}); org_sci_name = org_data.get("scientificName", "N/A"); org_common_name = org_data.get("commonName", "")
        org_display = f"{org_sci_name}" + (f" ({org_common_name})" if org_common_name else ""); seq_info = data.get("sequence", {}); seq_val = seq_info.get("value", "N/A"); length = seq_info.get("length", 0)
        status_val = data.get("entryAudit", {}).get("entryType", "N/A").replace("UniProtKB ", ""); existence_val = data.get("proteinExistence", "N/A").replace(": Evidence at ", ": ")
        score_val = data.get("annotationScore", "N/A"); mw_str = "N/A"
        if seq_val != "N/A" and length > 0:
            try:
                clean_seq = "".join(filter(lambda x: x in AMINO_ACID_NAMES, seq_val.upper()))
                if clean_seq: mw_str = f"{molecular_weight(clean_seq, seq_type='protein'):.2f} Da"
                else: mw_str = "Invalid sequence for MW"
            except: mw_str = "Error in MW calc"
        
        # --- Corrected Indentation ---
        comments_data = data.get("comments", [])
        func_comment = "N/A" 
        for c_item in comments_data:           
            if c_item.get("commentType") == "FUNCTION":
                texts = c_item.get("texts", []) 
                if texts:                   
                    func_comment = texts[0].get("value", "N/A") 
                    break                     
        # --- End Correction ---

        overview_md = (f"## {id_display} ({acc})\n[{acc} on UniProt]({link})\n\n**Protein:** {name}\n**Gene:** {gene_str}\n**Status:** {status_val}\n"
                       f"**Organism:** {org_display}\n**Length:** {length} aa\n**Existence:** {existence_val}\n**Score:** {score_val}/5\n**Calc. MW:** {mw_str}\n\n"
                       f"**Function Snippet:**\n{func_comment}\n\n**Sequence (first 100 aa):**\n`{seq_val[:100]}{'...' if len(seq_val) > 100 else ''}`\n\n--- \n*More details in other tabs.*")
        
        # Extract other info
        interactions_md = extract_interactions(data); pathways_md = extract_pathways(data)
        disease_md = extract_disease_info(data); publications_md = extract_publications(data)
        xref_md = extract_cross_references(data) 

        # Process plots
        aa_freq, aa_err = get_amino_acid_frequencies(seq_val); aa_plot_upd = empty_plot
        if aa_err: overview_md += f"\n\n**AA Freq Error:** {aa_err}" # Append error if needed
        elif aa_freq: 
            img_aa = plot_amino_acid_frequencies(aa_freq)
            if img_aa: aa_plot_upd = gr.update(value=img_aa, visible=True)
        
        seq_feat = extract_sequence_features(data); feat_plot_upd = empty_plot; feat_msg = ""
        if seq_feat and length > 0:
            img_feat = plot_sequence_features(length, seq_feat)
            if img_feat: feat_plot_upd = gr.update(value=img_feat, visible=True)
            else: feat_msg = "Could not generate feature plot."
        elif not seq_feat and length > 0 : feat_msg = "No relevant features found for plotting."
        
        # Return all outputs in the correct order for the UI
        return (overview_md, aa_plot_upd, feat_plot_upd, feat_msg, 
                pathways_md, interactions_md, disease_md, publications_md, xref_md)
    
    # Error Handling for the main function
    except requests.exceptions.HTTPError as e: 
        err_msg_http = f"Error: ID '{uniprot_id}' not found." if e.response.status_code == 404 else f"HTTP error: {e}"
        return (err_msg_http,) + outputs_on_error[1:]
    except Exception as e:
        return (f"Error: {str(e)[:150]}",) + outputs_on_error[1:]

# --- Gradio UI Definition ---
with gr.Blocks(theme=gr.themes.Glass()) as iface:
    gr.Markdown("# Protein Profile Viewer (v1.6.1 - Indent Fix)")
    gr.Markdown("Enter a UniProt ID directly, **OR** search for a protein/gene name below to find and copy its ID.")
    
    with gr.Group():
        gr.Markdown("### Find UniProt ID by Name/Keyword")
        search_term_input = gr.Textbox(label="Search Term (min 3 chars, type and wait)", placeholder="e.g., insulin, EGFR, P53")
        search_status_output = gr.Markdown() 
        search_results_output = gr.Dataset(
            label="Search Results (Select a row to copy Accession)", 
            headers=["Accession", "UniProtKB ID", "Protein Name", "Organism"],
            samples=[], samples_per_page=5 
        )
        selected_id_to_copy = gr.Textbox(label="Copy this Accession ID:", interactive=True, show_copy_button=True)

    gr.Markdown("---") 

    gr.Markdown("### View Protein Profile")
    with gr.Row():
        protein_id_input = gr.Textbox(label="Enter UniProt Accession ID", placeholder="Paste ID here or enter directly", scale=3)
        submit_button = gr.Button("Get Profile", scale=1, variant="primary")

    with gr.Tabs():
        with gr.TabItem("Overview"):
             gr.Markdown("### Protein Overview\nKey info: UniProt ID, name, gene, organism, function, sequence snippet & UniProt link.")
             overview_output = gr.Markdown()
        with gr.TabItem("Analysis Plots"):
             gr.Markdown("### Sequence Analysis Visualizations\nAmino acid composition and annotated features.")
             with gr.Column(): 
                 aa_freq_plot_output = gr.Image(label="Amino Acid Frequency Plot", type="pil", show_label=True, visible=False)
                 seq_features_plot_output = gr.Image(label="Sequence Features Plot", type="pil", show_label=True, visible=False)
                 seq_features_message_output = gr.Markdown() 
        with gr.TabItem("Functional Context"):
             gr.Markdown("### Pathways, Interactions & Disease\nBiological context.")
             with gr.Accordion("Biological Pathways", open=False): pathways_output = gr.Markdown()
             with gr.Accordion("Protein Interactions", open=False): interactions_output = gr.Markdown()
             with gr.Accordion("Disease Associations", open=False): disease_output = gr.Markdown()
        with gr.TabItem("Publications"):
             gr.Markdown("### Relevant Publications\nAssociated scientific literature.")
             publications_output = gr.Markdown()
        with gr.TabItem("Cross-references"):
             gr.Markdown("### Database Links\nLinks to other relevant databases.")
             xref_output = gr.Markdown()

    # --- Event Handlers ---
    search_term_input.change(fn=search_uniprot_by_name, inputs=search_term_input, outputs=[search_status_output, search_results_output])
    search_results_output.select(fn=update_copy_box, inputs=None, outputs=selected_id_to_copy, _js="(evt) => evt.value[0]") # Correct JS for Dataset select
    submit_button.click(
        fn=get_protein_info, inputs=protein_id_input,
        outputs=[overview_output, aa_freq_plot_output, seq_features_plot_output, seq_features_message_output, 
                 pathways_output, interactions_output, disease_output, publications_output, xref_output] 
    )
    gr.Examples(examples=[["P05067"], ["P00533"], ["Q9BYF1"], ["P0DP23"], ["P04637"]], inputs=protein_id_input)

if __name__ == "__main__":
    iface.launch()
