import gradio as gr
import requests # برای ارسال درخواست به API
from Bio.SeqUtils import molecular_weight # برای محاسبه وزن مولکولی

# آدرس پایه UniProt API برای دریافت اطلاعات پروتئین به فرمت JSON
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"

def get_protein_info(uniprot_id):
    """
    اطلاعات یک پروتئین را از UniProt API دریافت و پردازش می‌کند.
    """
    if not uniprot_id:
        return " لطفاً یک شناسه UniProt وارد کنید.", "", "", "", "", ""

    # ساخت URL کامل برای درخواست API
    url = UNIPROT_API_URL.format(accession=uniprot_id.strip().upper())

    try:
        response = requests.get(url)
        response.raise_for_status()  # اگر خطایی در درخواست بود (مثل 404)، اینجا متوقف می‌شود
        data = response.json()

        # استخراج اطلاعات اصلی
        primary_accession = data.get("primaryAccession", "N/A")
        protein_data = data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "N/A")
        if protein_data == "N/A" and data.get("proteinDescription", {}).get("submissionNames"):
            protein_data = data["proteinDescription"]["submissionNames"][0].get("fullName", {}).get("value", "N/A")

        gene_names_data = data.get("genes")
        gene_name_str = "N/A"
        if gene_names_data:
            gene_name_str = ", ".join([g.get("geneName", {}).get("value", "") for g in gene_names_data if g.get("geneName")])
            if not gene_name_str: # اگر فقط ORFName یا OrderedLocusName بود
                gene_name_str = ", ".join([g.get("orfNames", [{}])[0].get("value", "") for g in gene_names_data if g.get("orfNames")])


        organism_data = data.get("organism", {})
        organism_name = organism_data.get("scientificName", "N/A")
        if organism_data.get("commonName"):
            organism_name += f" ({organism_data.get('commonName')})"

        sequence_info = data.get("sequence", {})
        sequence = sequence_info.get("value", "N/A")
        length = sequence_info.get("length", 0)
        # mol_weight = sequence_info.get("molWeight", 0) # UniProt گاهی اینو مستقیم میده

        # محاسبه وزن مولکولی با BioPython اگر در API نبود یا برای دقت بیشتر
        # BioPython انتظار داره توالی فقط شامل حروف استاندارد آمینواسید باشه
        # پس باید حروف غیر استاندارد (مثل X, B, Z, U, O) رو مدیریت کنیم اگر وجود دارن
        # برای سادگی فعلا فرض میکنیم توالی استاندارد است یا خطا میدهیم
        calculated_mol_weight_str = "N/A"
        if sequence != "N/A":
            try:
                # حذف کاراکترهای غیر استاندارد برای محاسبه (ساده سازی)
                cleaned_sequence = "".join(filter(lambda x: x in "ACDEFGHIKLMNPQRSTVWY", sequence.upper()))
                if cleaned_sequence:
                     calculated_mol_weight = molecular_weight(cleaned_sequence, seq_type="protein")
                     calculated_mol_weight_str = f"{calculated_mol_weight:.2f} Da"
                else:
                    calculated_mol_weight_str = "توالی نامعتبر برای محاسبه"
            except Exception as e:
                calculated_mol_weight_str = f"خطا در محاسبه: {e}"


        # اطلاعات عملکرد (بخش اول توضیحات عملکرد)
        comments = data.get("comments", [])
        function_comment = "N/A"
        for comment in comments:
            if comment.get("commentType") == "FUNCTION":
                function_texts = comment.get("texts", [])
                if function_texts:
                    function_comment = function_texts[0].get("value", "N/A")
                    break # فقط اولین پاراگراف عملکرد

        return (
            f"**شناسه اصلی:** {primary_accession}\n"
            f"**نام پروتئین:** {protein_data}\n"
            f"**نام ژن:** {gene_name_str}\n"
            f"**ارگانیسم:** {organism_name}\n"
            f"**طول توالی:** {length} آمینو اسید\n"
            f"**وزن مولکولی (محاسبه شده):** {calculated_mol_weight_str}\n\n"
            f"**بخشی از عملکرد:**\n{function_comment}\n\n"
            f"**توالی آمینواسیدی (100 حرف اول):**\n`{sequence[:100]}{'...' if len(sequence) > 100 else ''}`"
        )

    except requests.exceptions.HTTPError as http_err:
        if response.status_code == 404:
            return f"خطا: پروتئین با شناسه '{uniprot_id}' یافت نشد. لطفاً شناسه را بررسی کنید.", "", "", "", "", ""
        return f"خطای HTTP رخ داد: {http_err}", "", "", "", "", ""
    except requests.exceptions.RequestException as req_err:
        return f"خطای شبکه رخ داد: {req_err}", "", "", "", "", ""
    except Exception as e:
        return f"یک خطای غیرمنتظره رخ داد: {e}\n بررسی کنید شناسه معتبر باشد.", "", "", "", "", ""

# تعریف رابط کاربری Gradio
iface = gr.Interface(
    fn=get_protein_info,
    inputs=gr.Textbox(label="شناسه UniProt را وارد کنید (مثال: P05067 یا INS_HUMAN)", placeholder="مثلاً P0DP23"),
    outputs=gr.Markdown(label="اطلاعات پروتئین"), # استفاده از Markdown برای نمایش بهتر متن فرمت شده
    title="نمایشگر پروفایل پروتئین (نسخه پایه)",
    description="شناسه یک پروتئین از UniProt را وارد کنید تا اطلاعات اولیه آن نمایش داده شود.",
    examples=[["P05067"], ["INS_HUMAN"], ["CYC_HUMAN"], ["P12345"]], # مثال‌هایی برای تست سریع
    allow_flagging='never' # برای غیرفعال کردن دکمه Flag
)

if __name__ == "__main__":
    iface.launch() # برای اجرای محلی
    # برای Hugging Face Spaces، فقط iface.launch() کافی است و نیازی به بلاک __main__ نیست
    # اما برای اجرای محلی خوب است که این بلاک وجود داشته باشد.
