import time
import os
from dotenv import load_dotenv
from google import genai
from Bio import Entrez
from fpdf import FPDF

# 1. Carrega as senhas do arquivo .env (Seguranca)
load_dotenv()

# --- CONFIGURAÇÕES ---
MINHA_CHAVE = os.getenv("GEMINI_API_KEY")
Entrez.email = os.getenv("ENTREZ_EMAIL")

if not MINHA_CHAVE:
    print("\n[ERRO] A chave GEMINI_API_KEY nao foi encontrada no arquivo .env!")
    print("Verifique se o arquivo .env existe e se o nome da chave esta correto.")
    exit()

client = genai.Client(api_key=MINHA_CHAVE)

# --- BANCO DE DADOS DE CLIENTES (SaaS Model) ---
clientes = {
    "1": {"nome": "Dra. Ana Paula", "clinica": "NeuroVida", "especialidade": "Psiquiatria", "plano": "Premium", "limite": 5},
    "2": {"nome": "Dr. Carlos Alberto", "clinica": "Clinic-Ar", "especialidade": "Cardiologia", "plano": "Basico", "limite": 1},
    "3": {"nome": "Cunhada Querida", "clinica": "Consultorio Particular", "especialidade": "Psiquiatria", "plano": "Premium", "limite": 3}
}

class PDF_Personalizado(FPDF):
    def __init__(self, cliente_info):
        super().__init__()
        self.cliente = cliente_info

    def header(self):
        self.set_fill_color(0, 51, 102) 
        self.rect(0, 0, 210, 40, 'F')
        self.set_font("helvetica", 'B', 20)
        self.set_text_color(255, 255, 255)
        self.cell(0, 15, "MEDICAL IN-SIGHT PREMIUM", new_x="LMARGIN", new_y="NEXT", align='C')
        self.set_font("helvetica", 'I', 10)
        self.cell(0, 5, f"Relatorio Exclusivo: {self.cliente['clinica']}", new_x="LMARGIN", new_y="NEXT", align='C')
        self.ln(20)

    def footer(self):
        self.set_y(-15)
        self.set_font("helvetica", 'I', 8)
        self.set_text_color(150, 150, 150)
        self.cell(0, 10, f"Preparado para {self.cliente['nome']} | Pagina {self.page_no()}", align='C')

def buscar_pubmed(tema, limite):
    print(f"\n[1/3] Buscando {limite} evidencias cientificas sobre: {tema}...")
    handle = Entrez.esearch(db="pubmed", term=tema, retmax=limite, sort="date")
    record = Entrez.read(handle)
    handle.close()
    
    artigos = []
    for id_artigo in record.get("IdList", []):
        fetch = Entrez.efetch(db="pubmed", id=id_artigo, rettype="abstract", retmode="text")
        artigos.append({
            "texto": fetch.read(),
            "link": f"https://pubmed.ncbi.nlm.nih.gov/{id_artigo}/"
        })
        fetch.close()
    return artigos

def main():
    print("\n" + "="*40)
    print("   BEM-VINDO AO PORTAL MEDICAL IN-SIGHT   ")
    print("="*40)
    
    for id_c, info in clientes.items():
        print(f"{id_c}. {info['nome']} ({info['clinica']})")
    
    escolha = input("\nSelecione o Cliente para gerar o boletim: ")
    
    if escolha in clientes:
        user = clientes[escolha]
        print(f"\nIniciando curadoria para {user['nome']}...")
        
        artigos_brutos = buscar_pubmed(user['especialidade'], user['limite'])
        
        if artigos_brutos:
            print("[2/3] Gemini 2.0 analisando e vinculando fontes...")
            texto_para_ia = "\n\n".join([a['texto'] for a in artigos_brutos])
            
            prompt = f"""
            Aja como um curador cientifico para o(a) {user['nome']}.
            Traduza e resuma os artigos abaixo.
            O tom deve ser profissional e focado na especialidade {user['especialidade']}.
            IMPORTANTE: Termine cada resumo individual com a palavra exata: [FONTE]
            
            ARTIGOS:
            {texto_para_ia[:8000]}
            """
            
            response = client.models.generate_content(model="gemini-2.0-flash", contents=prompt)
            
            print("[3/3] Renderizando PDF com Links Centralizados...")
            pdf = PDF_Personalizado(user)
            pdf.add_page()
            
            pdf.set_font("helvetica", 'B', 12)
            pdf.set_text_color(0, 51, 102)
            pdf.cell(0, 10, f"DATA: {time.strftime('%d/%m/%Y')} | DESTINATARIO: {user['nome'].upper()}", new_x="LMARGIN", new_y="NEXT")
            pdf.ln(5)
            
            partes_do_texto = response.text.split("[FONTE]")
            
            for i, parte in enumerate(partes_do_texto):
                if parte.strip():
                    # Escreve o resumo
                    pdf.set_font("helvetica", size=11)
                    pdf.set_text_color(0, 0, 0)
                    texto_seguro = parte.encode('latin-1', 'replace').decode('latin-1')
                    pdf.multi_cell(0, 8, txt=texto_seguro)
                    
                    # Escreve o link centralizado (DENTRO do IF para alinhar com o texto)
                    if i < len(artigos_brutos):
                        pdf.ln(2)
                        pdf.set_font("helvetica", 'B', 10)
                        pdf.set_text_color(0, 0, 255)
                        link_url = artigos_brutos[i]['link']
                        pdf.cell(0, 10, txt="--- CLIQUE AQUI PARA LER O ESTUDO COMPLETO ---", 
                                 link=link_url, new_x="LMARGIN", new_y="NEXT", align='C')
                        pdf.ln(10) # Espaço entre artigos
            
            arquivo = f"Boletim_{user['nome'].replace(' ', '_')}.pdf"
            pdf.output(arquivo)
            print(f"\n--- SUCESSO! PDF gerado: {arquivo} ---")
        else:
            print("Nenhum artigo novo encontrado.")
    else:
        print("Opcao invalida.")

if __name__ == "__main__":
    main()