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
import pandas as pd # For potential use with Dataset if needed, though list of lists works

UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

AMINO_ACID_NAMES = { # ... (same as before) ... 
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic acid',
    'C': 'Cysteine', 'Q': 'Glutamine', 'E': 'Glutamic acid', 'G': 'Glycine',
    'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
    'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
    'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine'
}
STANDARD_AMINO_ACIDS_ORDER = "ARNDCQEGHILKMFPSTWYV" 

# --- Helper Functions (Plotting, Feature/Interaction/etc. Extraction - keep as is) ---
def get_amino_acid_frequencies(sequence): # ... (same as before) ...
    if not sequence or sequence == "N/A": return None, "Sequence not available for analysis."
    cleaned_sequence = "".join(filter(lambda x: x in AMINO_ACID_NAMES, sequence.upper()))
    if not cleaned_sequence: return None, "No valid amino acids found for counting."
    counts = Counter(cleaned_sequence)
    frequencies = {aa: counts.get(aa, 0) for aa in STANDARD_AMINO_ACIDS_ORDER}
    return frequencies, None

def plot_amino_acid_frequencies(frequencies): # ... (same as before) ...
    if not frequencies: return None
    ordered_keys = [key for key in STANDARD_AMINO_ACIDS_ORDER if key in frequencies]
    labels = [f"{aa}: {AMINO_ACID_NAMES.get(aa, aa)}" for aa in ordered_keys]
    values = [frequencies[aa] for aa in ordered_keys]
    fig, ax = plt.subplots(figsize=(12, 7)); ax.bar(labels, values, color='skyblue')
    ax.set_xlabel("Amino Acid"); ax.set_ylabel("Frequency"); ax.set_title("AA Freq Plot")
    plt.xticks(rotation=75, ha="right", fontsize=8); plt.tight_layout()
    buf = io.BytesIO(); plt.savefig(buf, format='png'); buf.seek(0)
    img = Image.open(buf); plt.close(fig); return img

def extract_sequence_features(uniprot_data): # ... (same as before) ...
    features_of_interest_uppercase = {"DOMAIN": "blue", "MOTIF": "green", "ACTIVE_SITE": "red", "BINDING_SITE": "orange", "MOD_RES": "purple", "HELIX": "cyan", "STRAND": "magenta", "TURN": "gold"}
    extracted = []
    if "features" in uniprot_data and uniprot_data["features"]:
        for item in uniprot_data["features"]:
            type_raw = item.get("type"); loc_obj = item.get("location", {})
            if not isinstance(type_raw, str): continue
            type_norm = type_raw.strip().upper()
            if type_norm in features_of_interest_uppercase:
                try:
                    b_str, e_str = None, None; s_node, e_node, p_node = loc_obj.get("start"), loc_obj.get("end"), loc_obj.get("position")
                    if s_node and isinstance(s_node, dict) and "value" in s_node: b_str = str(s_node["value"])
                    if e_node and isinstance(e_node, dict) and "value" in e_node: e_str = str(e_node["value"])
                    if p_node and isinstance(p_node, dict) and "value" in p_node:
                        p_str = str(p_node["value"]);
                        if b_str is None: b_str = p_str
                        if e_str is None: e_str = p_str
                        if "start" not in loc_obj and "end" not in loc_obj: b_str, e_str = p_str, p_str
                    if b_str is None or e_str is None: continue
                    b_pos, e_pos = int(b_str), int(e_str)
                    if b_pos > e_pos: continue
                    extracted.append({"type": type_raw, "begin": b_pos, "end": e_pos, "description": item.get("description", type_raw), "color": features_of_interest_uppercase[type_norm]})
                except: continue
    return extracted

def plot_sequence_features(sequence_length, features): # ... (same as before) ...
    if not features or sequence_length == 0: return None
    fig, ax = plt.subplots(figsize=(12, max(3, len(features) * 0.4) + 1.5)) 
    ax.set_xlim(0, sequence_length); ax.set_xlabel("AA Position"); ax.set_yticks([]) 
    ax.set_title("Sequence Features"); leg_h = {} 
    y_pos, leg_set, b_h = 0, set(), 0.8 
    for feat in sorted(features, key=lambda x: x["begin"]):
        b, e, c = feat["begin"], feat["end"], feat["color"]; w = max(1, e - b + 1) 
        ax.barh(y_pos, w, height=b_h, left=b -1, color=c, edgecolor='black', alpha=0.7)
        if feat['type'] not in leg_set: leg_h[feat['type']] = plt.Rectangle((0, 0), 1, 1, fc=c, alpha=0.7); leg_set.add(feat['type'])
        y_pos += 1 
    if y_pos > 0: ax.set_ylim(-0.5, y_pos -1 + b_h/2 + 0.5) 
    else: plt.close(fig); return None
    if leg_h: ax.legend(leg_h.values(), leg_h.keys(), title="Feature Types", bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    plt.tight_layout(rect=[0, 0, 0.83, 0.96])
    buf = io.BytesIO(); plt.savefig(buf, format='png'); buf.seek(0); img = Image.open(buf); plt.close(fig); return img

def extract_interactions(uniprot_data): # ... (same as before) ...
    interactions = []
    if "comments" in uniprot_data:
        for c in uniprot_data["comments"]:
            if c.get("commentType") == "INTERACTION" and "interactions" in c:
                for i_entry in c["interactions"]:
                    i1_acc = i_entry.get("interactantOne", {}).get("uniProtKBAccession")
                    i2_acc = i_entry.get("interactantTwo", {}).get("uniProtKBAccession")
                    i2_gene = i_entry.get("interactantTwo", {}).get("geneName")
                    if i1_acc == uniprot_data.get("primaryAccession") and i2_acc:
                        partner = f"{i2_gene} ({i2_acc})" if i2_gene else i2_acc
                        interactions.append(f"- Interacts with: **{partner}**")
    return "\n".join(sorted(interactions)) if interactions else "No interaction partners listed in comments."

def extract_pathways(uniprot_data): # ... (same as before) ...
    pathways = []; dbs = { "KEGG": "...", "Reactome": "..." } # URLs omitted for brevity
    if "uniProtKBCrossReferences" in uniprot_data:
        for xref in uniprot_data["uniProtKBCrossReferences"]:
            db = xref.get("database"); pid = xref.get("id")
            if db in dbs and pid:
                desc = pid # Default
                if "properties" in xref:
                    for p in xref["properties"]:
                        if p.get("key") in ["PathwayName", "Description"]: desc = f"{p.get('value')} ({pid})"; break
                link = dbs[db] + pid
                pathways.append(f"- [{desc}]({link}) ({db})")
    return "\n".join(sorted(list(set(pathways)))) if pathways else "No KEGG/Reactome pathway info."

def extract_disease_info(uniprot_data): # ... (same as before) ...
    diseases = []
    if "comments" in uniprot_data:
        for c in uniprot_data["comments"]:
            if c.get("commentType") == "DISEASE" and "disease" in c:
                d_entry = c["disease"]; d_name = d_entry.get("diseaseId", "?"); desc = d_entry.get("description", "N/A")
                mim = None;
                if "diseaseCrossReference" in d_entry and d_entry["diseaseCrossReference"].get("database") == "MIM": mim = d_entry["diseaseCrossReference"].get("id")
                d_md = f"**{d_name}**" + (f" (MIM: [{mim}](https://omim.org/entry/{mim}))" if mim else "")
                d_md += f"\n   - *Desc:* {desc}\n"; note_val = c.get("note", {}).get("texts", [{}])[0].get("value")
                if note_val: d_md += f"   - *Note:* {note_val}\n"
                diseases.append(d_md)
    return "\n---\n".join(sorted(diseases)) if diseases else "No disease association info."

def extract_publications(uniprot_data): # ... (same as before) ...
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
            md = f"**{i + 1}. {title}**\n   - *{authors}*\n   - *{j}" + (f", {v}" if v else "") + (f":{f}" if f else "") + (f"-{l}" if l else "") + (f" ({d})" if d else "") + "*\n"
            if pmid: md += f"   - [PubMed {pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)\n"
            if doi: md += f"   - [DOI {doi}](https://doi.org/{doi})\n"
            pubs.append(md)
    return "\n---\n".join(pubs) if pubs else "No publication info."

def extract_cross_references(uniprot_data): # ... (same as before) ...
    xrefs = []; dbs = {"Ensembl": "...", "GeneID": "...", "RefSeq": "...", "GO": "...", "InterPro": "...", "Pfam": "...", "PDB": "...", "KEGG": "...", "Reactome": "..."} # URLs omitted
    grouped = {db: [] for db in dbs};
    if "uniProtKBCrossReferences" in uniprot_data:
        for xr in uniprot_data["uniProtKBCrossReferences"]:
            db = xr.get("database"); xid = xr.get("id")
            if db in dbs and xid:
                url, txt = None, xid
                # Simplified logic for brevity
                if db == "GO": txt = xr.get("properties", [{}])[0].get("value", xid) + f" ({xid})"
                url = dbs[db] + xid # Simplified URL creation
                if url: grouped[db].append(f"[{txt}]({url})")
        for db, links in grouped.items():
            if links: xrefs.append(f"**{db}:** " + ", ".join(sorted(list(set(links)))))
    return "\n".join(xrefs) if xrefs else "No selected cross-references."


# --- NEW Search Function ---
def search_uniprot_by_name(search_term, result_limit=5):
    if not search_term or len(search_term) < 3: 
        return "Enter at least 3 characters.", [] # Return empty list for dataset
    params = { "query": f'({search_term}) AND (reviewed:true)', "fields": "accession,id,protein_name,organism_name", "format": "json", "size": result_limit }
    status_md = f"### Search Results for '{search_term}':\n"
    dataset_data = []
    try:
        response = requests.get(UNIPROT_SEARCH_URL, params=params); response.raise_for_status(); data = response.json()
        results = data.get("results")
        if not results: status_md += "No reviewed entries found."
        else:
            status_md += f"*Found {len(results)}. Select a row below, then copy the ID from the 'Copy this ID' box.*\n---"
            for entry in results:
                acc = entry.get("primaryAccession", "N/A"); uid = entry.get("uniProtkbId", ""); name = "N/A"
                if entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"):
                    name = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
                elif entry.get("proteinDescription", {}).get("submissionNames"): name = entry["proteinDescription"]["submissionNames"][0].get("fullName",{}).get("value", "N/A")
                org = entry.get("organism", {}).get("scientificName", "N/A")
                dataset_data.append([acc, uid, name, org])
    except Exception as e: status_md += f"\n**Search Error:** {e}"
    return status_md, dataset_data

# --- NEW Function to update copy box ---
def update_copy_box(evt: gr.SelectData):
    # evt.value should be the value of the first cell in the selected row (Accession ID)
    if evt.value:
        return evt.value
    return "" # Return empty string if no value (e.g., selection cleared)

# --- Main Function to Get All Info ---
def get_protein_info(uniprot_id):
    # ... (Implementation from v1.5.1, returns 9 outputs) ...
    empty_plot = gr.update(value=None, visible=False); empty_str = ""; err_msg = "Error"
    outputs_on_error = (err_msg, empty_plot, empty_plot, empty_str, empty_str, empty_str, empty_str, empty_str, empty_str) 
    if not uniprot_id: return ("Enter UniProt ID.",) + outputs_on_error[1:]
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
        comments_data = data.get("comments", []); func_comment = "N/A" 
        for c_item in comments_data:
            if c_item.get("commentType") == "FUNCTION": texts = c_item.get("texts", []);
                if texts: func_comment = texts[0].get("value", "N/A"); break
        overview_md = (f"## {id_display} ({acc})\n[{acc} on UniProt]({link})\n\n**Protein:** {name}\n**Gene:** {gene_str}\n**Status:** {status_val}\n"
                       f"**Organism:** {org_display}\n**Length:** {length} aa\n**Existence:** {existence_val}\n**Score:** {score_val}/5\n**Calc. MW:** {mw_str}\n\n"
                       f"**Function Snippet:**\n{func_comment}\n\n**Sequence (first 100 aa):**\n`{seq_val[:100]}{'...' if len(seq_val) > 100 else ''}`\n\n--- \n*More details in other tabs.*")
        interactions_md = extract_interactions(data); pathways_md = extract_pathways(data)
        disease_md = extract_disease_info(data); publications_md = extract_publications(data); xref_md = extract_cross_references(data) 
        aa_freq, aa_err = get_amino_acid_frequencies(seq_val); aa_plot_upd = empty_plot
        if aa_err: overview_md += f"\n\n**AA Freq Error:** {aa_err}"
        elif aa_freq: img_aa = plot_amino_acid_frequencies(aa_freq);
            if img_aa: aa_plot_upd = gr.update(value=img_aa, visible=True)
        seq_feat = extract_sequence_features(data); feat_plot_upd = empty_plot; feat_msg = ""
        if seq_feat and length > 0:
            img_feat = plot_sequence_features(length, seq_feat)
            if img_feat: feat_plot_upd = gr.update(value=img_feat, visible=True)
            else: feat_msg = "Could not generate feature plot."
        elif not seq_feat and length > 0 : feat_msg = "No relevant features found for plotting."
        return (overview_md, aa_plot_upd, feat_plot_upd, feat_msg, pathways_md, interactions_md, disease_md, publications_md, xref_md)
    except requests.exceptions.HTTPError as e: 
        err_msg_http = f"Error: ID '{uniprot_id}' not found." if e.response.status_code == 404 else f"HTTP error: {e}"
        return (err_msg_http,) + outputs_on_error[1:]
    except Exception as e: return (f"Error: {str(e)[:150]}",) + outputs_on_error[1:]

# --- Gradio UI Definition ---
with gr.Blocks(theme=gr.themes.Glass()) as iface:
    gr.Markdown("# Protein Profile Viewer (v1.6 - Live Search & Copy)")
    gr.Markdown("Enter a UniProt ID directly, **OR** search for a protein/gene name below to find and copy its ID.")
    
    with gr.Group():
        gr.Markdown("### Find UniProt ID by Name/Keyword")
        search_term_input = gr.Textbox(label="Search Term (min 3 chars, type and wait)", placeholder="e.g., insulin, EGFR, P53")
        search_status_output = gr.Markdown() # To show "Searching..." or "No results"
        search_results_output = gr.Dataset(
            label="Search Results (Select a row to copy Accession)", 
            headers=["Accession", "UniProtKB ID", "Protein Name", "Organism"],
            samples=[], # Start empty
            samples_per_page=5 
        )
        selected_id_to_copy = gr.Textbox(
            label="Copy this Accession ID:", 
            interactive=True, # Allows user interaction (like copying)
            show_copy_button=True # Explicitly show the copy button
        )

    gr.Markdown("---") 

    gr.Markdown("### View Protein Profile")
    with gr.Row():
        # Note: User now copies from selected_id_to_copy and pastes here
        protein_id_input = gr.Textbox(label="Enter UniProt Accession ID", placeholder="Paste ID here from search results or enter directly", scale=3)
        submit_button = gr.Button("Get Profile", scale=1, variant="primary")

    with gr.Tabs():
        with gr.TabItem("Overview"):
            gr.Markdown("### Protein Overview\nKey information...")
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
    # Live search trigger (when text changes in search box)
    search_term_input.change(
        fn=search_uniprot_by_name,
        inputs=search_term_input,
        outputs=[search_status_output, search_results_output] 
        # Note: live=True is implicit with .change() for Textbox by default
        # We could add debounce here if Gradio supports it in the future or use a more complex setup.
    )

    # When a row is selected in the Dataset, update the copy box
    search_results_output.select(
        fn=update_copy_box,
        inputs=None, # Input evt is handled implicitly by SelectData
        outputs=selected_id_to_copy,
        _js="(evt) => evt.value" # Pass the value of the first cell (Accession) directly
    )

    # Main profile submit button action
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
