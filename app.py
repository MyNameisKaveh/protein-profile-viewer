import gradio as gr
import requests
from Bio.SeqUtils import molecular_weight
from collections import Counter
import matplotlib.pyplot as plt
import io
# import base64 # دیگه برای این روش لازم نیست
from PIL import Image # برای کار با آبجکت تصویر Pillow

UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"

def get_amino_acid_frequencies(sequence):
    # ... (این تابع بدون تغییر باقی می‌مونه) ...
    if not sequence or sequence == "N/A":
        return None, "توالی برای تحلیل موجود نیست."
    standard_amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    cleaned_sequence = "".join(filter(lambda x: x in standard_amino_acids, sequence.upper()))
    if not cleaned_sequence:
        return None, "توالی معتبر برای شمارش آمینواسیدها یافت نشد."
    counts = Counter(cleaned_sequence)
    frequencies = {aa: counts.get(aa, 0) for aa in standard_amino_acids}
    return frequencies, None

def plot_amino_acid_frequencies(frequencies):
    """
    نمودار میله‌ای فراوانی آمینواسیدها را رسم می‌کند و به صورت آبجکت تصویر Pillow برمی‌گرداند.
    """
    if not frequencies:
        return None

    names = list(frequencies.keys())
    values = list(frequencies.values())

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(names, values, color='skyblue')
    ax.set_xlabel("آمینواسید")
    ax.set_ylabel("فراوانی")
    ax.set_title("نمودار فراوانی آمینواسیدها")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    img = Image.open(buf) # تبدیل بافر به آبجکت تصویر Pillow
    plt.close(fig) # مهم: بستن شکل برای آزاد کردن حافظه
    # buf.close() # بافر رو هم می‌بندیم (اختیاری، Image.open ممکنه خودش این کار رو بکنه)
    return img


def get_protein_info(uniprot_id):
    if not uniprot_id:
        return " لطفاً یک شناسه UniProt وارد کنید.", None

    url = UNIPROT_API_URL.format(accession=uniprot_id.strip().upper())

    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()

        # ... (بخش استخراج اطلاعات متنی بدون تغییر باقی می‌مونه) ...
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
                    calculated_mol_weight_str = "توالی نامعتبر"
            except Exception: 
                calculated_mol_weight_str = "خطا در محاسبه"

        comments = data.get("comments", [])
        function_comment = "N/A"
        for comment in comments:
            if comment.get("commentType") == "FUNCTION":
                function_texts = comment.get("texts", [])
                if function_texts:
                    function_comment = function_texts[0].get("value", "N/A")
                    break
        
        main_info_md = (
            f"**شناسه اصلی:** {primary_accession}\n"
            f"**نام پروتئین:** {protein_name}\n"
            f"**نام ژن:** {gene_name_str}\n"
            f"**ارگانیسم:** {organism_name}\n"
            f"**طول توالی:** {length} آمینو اسید\n"
            f"**وزن مولکولی (محاسبه شده):** {calculated_mol_weight_str}\n\n"
            f"**بخشی از عملکرد:**\n{function_comment}\n\n"
            f"**توالی آمینواسیدی (100 حرف اول):**\n`{sequence[:100]}{'...' if len(sequence) > 100 else ''}`"
        )
        
        aa_frequencies, freq_error = get_amino_acid_frequencies(sequence)
        aa_plot_pil_image = None # حالا یک آبجکت تصویر Pillow خواهد بود
        if freq_error:
            main_info_md += f"\n\n**تحلیل فراوانی آمینواسید:** {freq_error}"
        elif aa_frequencies:
            aa_plot_pil_image = plot_amino_acid_frequencies(aa_frequencies)
        
        return main_info_md, aa_plot_pil_image

    # ... (بخش مدیریت خطا بدون تغییر باقی می‌مونه) ...
    except requests.exceptions.HTTPError as http_err:
        status_code = http_err.response.status_code if http_err.response is not None else "N/A"
        if status_code == 404:
            return f"خطا: پروتئین با شناسه '{uniprot_id}' یافت نشد. لطفاً شناسه را بررسی کنید.", None
        return f"خطای HTTP رخ داد: {http_err} (کد وضعیت: {status_code})", None
    except requests.exceptions.RequestException as req_err:
        return f"خطای شبکه رخ داد: {req_err}", None
    except Exception as e:
        return f"یک خطای غیرمنتظره در پردازش داده‌ها رخ داد. لطفاً شناسه را بررسی کنید و دوباره تلاش نمایید.", None

# تعریف رابط کاربری Gradio
# outputs_list بدون تغییر، چون gr.Image(type="pil") آبجکت Pillow رو قبول می‌کنه
outputs_list = [
    gr.Markdown(label="اطلاعات پروتئین"),
    gr.Image(label="نمودار فراوانی آمینواسیدها", type="pil", show_label=True) 
]

# اصلاح پارامتر allow_flagging به flagging_options (این هشدار هم در لاگ بود)
iface = gr.Interface(
    fn=get_protein_info,
    inputs=gr.Textbox(label="شناسه UniProt را وارد کنید (مثال: P05067 یا INS_HUMAN)", placeholder="مثلاً P0DP23"),
    outputs=outputs_list,
    title="نمایشگر پروفایل پروتئین (نسخه بهبود یافته)",
    description="شناسه یک پروتئین از UniProt را وارد کنید تا اطلاعات اولیه و نمودار فراوانی آمینواسیدهای آن نمایش داده شود.",
    examples=[["P05067"], ["INS_HUMAN"], ["CYC_HUMAN"], ["P12345"]],
    flagging_options=None # جایگزین allow_flagging='never' برای نسخه‌های جدید Gradio
)

if __name__ == "__main__":
    iface.launch()
