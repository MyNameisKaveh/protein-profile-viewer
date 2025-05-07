---
title: Protein Profile Viewer
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: gradio
sdk_version: "5.29.0"
app_file: app.py
pinned: false
---
# Protein Profile Viewer

[![Hugging Face Spaces](https://img.shields.io/badge/%F0%9F%A4%97%20Hugging%20Face-Spaces-blue)](https://huggingface.co/spaces/Andolinism/protein-profile-viewer) 
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0) 

An interactive web application built with Python and Gradio to fetch and display comprehensive information about a protein using its UniProt ID.

## Features

This application retrieves data directly from the UniProt API ([rest.uniprot.org](https://rest.uniprot.org/)) and presents it in an organized, user-friendly interface with multiple tabs:

*   **Overview Tab:**
    *   Displays fundamental protein details: UniProt Accession ID (with link), Protein Name, Gene Name(s), Organism, Sequence Length, Calculated Molecular Weight, Protein Existence Evidence, Annotation Score.
    *   Shows a snippet of the protein's function summary from UniProt comments.
    *   Provides the first 100 amino acids of the protein sequence.
*   **Analysis Plots Tab:**
    *   **Amino Acid Frequency Plot:** A bar chart visualizing the frequency of each standard amino acid in the protein sequence. Labels include the full amino acid name (e.g., "A: Alanine").
    *   **Sequence Features Plot:** A graphical representation showing the location and type of annotated features (Domains, Motifs, Active Sites, Secondary Structures like Helices and Strands, etc.) along the protein sequence. Includes a color-coded legend.
*   **Functional Context Tab:** (Information organized within accordions)
    *   **Biological Pathways:** Lists pathways from KEGG and Reactome databases that the protein is involved in, with direct links to the respective pathway pages.
    *   **Protein Interactions:** Displays known protein interaction partners as listed in UniProt comments.
    *   **Disease Associations:** Shows diseases linked to the protein, including descriptions and links to OMIM where available.
*   **Publications Tab:**
    *   Lists relevant scientific publications associated with the UniProt entry, including title, authors, journal details, and links to PubMed and DOI (if available).
*   **Cross-references Tab:**
    *   Provides direct links to the protein's entry in other major biological databases, such as Ensembl, NCBI Gene, RefSeq, Gene Ontology (GO), InterPro, Pfam, and PDB.

## Technologies Used

*   **Backend:** Python 3
*   **Web Framework/UI:** Gradio (`gradio`)
*   **Data Fetching:** `requests` (for UniProt API)
*   **Bioinformatics Calculations:** `BioPython` (`Bio.SeqUtils.molecular_weight`)
*   **Plotting:** `matplotlib`
*   **Image Handling:** `Pillow` (PIL Fork)
*   **Deployment:** Hugging Face Spaces
*   **CI/CD:** GitHub Actions (for automatic syncing from GitHub to HF Spaces)

## Getting Started

### Live Demo

You can access the live application hosted on Hugging Face Spaces:
**[https://huggingface.co/spaces/Andolinism/protein-profile-viewer](https://huggingface.co/spaces/Andolinism/protein-profile-viewer)**

### Local Installation and Usage (Optional)

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/MyNameisKaveh/protein-profile-viewer.git 
    cd protein-profile-viewer
    ```

2.  **Create and activate a virtual environment (recommended):**
    ```bash
    python -m venv venv
    # On Windows:
    # venv\Scripts\activate
    # On macOS/Linux:
    # source venv/bin/activate
    ```

3.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

4.  **Run the application:**
    ```bash
    python app.py
    ```
    The application will typically be available at `http://127.0.0.1:7860` in your web browser.

## How It Works

1.  The user enters a UniProt Accession ID (e.g., `P00533`).
2.  The Gradio interface sends the ID to the backend Python function (`get_protein_info`).
3.  The backend function constructs the UniProt API URL and sends a GET request using the `requests` library.
4.  The JSON response from UniProt is parsed.
5.  Various helper functions (`extract_...`, `plot_...`) process the JSON data to extract specific information (overview details, interactions, pathways, features, etc.) and generate plots using `matplotlib`.
6.  The processed information and plots are formatted as Markdown strings or PIL Image objects.
7.  These outputs are returned to the Gradio interface, which updates the content of the corresponding tabs and components.

## Future Enhancements / Potential Improvements

*   **3D Structure Viewer:** Re-integrate an interactive 3D viewer (e.g., NGL Viewer, Mol*) to display PDB structures when available. (Currently removed due to cross-platform/mobile rendering issues).
*   **More Detailed Interactions:** Parse interaction comments more deeply or integrate with dedicated interaction databases (e.g., IntAct API) to show interaction types and evidence.
*   **Clickable Features:** Make elements in the Sequence Features Plot clickable to show more details about a specific feature.
*   **Sequence Alignment/BLAST:** Add functionality to perform basic sequence similarity searches.
*   **Advanced Filtering/Search:** Allow filtering of interactions, pathways, or features.
*   **Error Handling:** More granular error handling and user feedback for API issues or data parsing problems.
*   **Caching:** Implement caching for API responses to speed up repeated lookups.
*   **Mobile Responsiveness:** Further investigate CSS/JS solutions to improve tab display on mobile if possible within Gradio constraints.

## Contributing

Contributions are welcome! If you find any issues or have suggestions for improvements, please open an issue or submit a pull request.

1.  Fork the repository.
2.  Create a new branch (`git checkout -b feature/your-feature-name`).
3.  Make your changes.
4.  Commit your changes (`git commit -am 'Add some feature'`).
5.  Push to the branch (`git push origin feature/your-feature-name`).
6.  Open a Pull Request.

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.
<!-- Make sure you have a LICENSE file with the Apache 2.0 text -->

## Acknowledgements

*   [UniProt](https://www.uniprot.org/) for providing the comprehensive protein data API.
*   [Gradio](https://www.gradio.app/) for the easy-to-use Python web UI framework.
*   [Hugging Face](https://huggingface.co/) for the free Spaces hosting platform.
*   The developers of BioPython, Matplotlib, Requests, and Pillow.
