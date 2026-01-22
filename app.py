import time
import os
import smtplib
from email.message import EmailMessage
from dotenv import load_dotenv
from google import genai
from Bio import Entrez
from fpdf import FPDF

# 1. Carrega as configurações de segurança do arquivo .env
load_dotenv()

# --- CONFIGURAÇÕES ---
MINHA_CHAVE = os.getenv("GEMINI_API_KEY")
Entrez.email = os.getenv("ENTREZ_EMAIL")
EMAIL_DE = os.getenv("EMAIL_REMETENTE")
SENHA_DE = os.getenv("EMAIL_SENHA")

if not MINHA_CHAVE or not EMAIL_DE or not SENHA_DE:
    print("\n[ERRO] Verifique se as chaves no arquivo .env estão preenchidas corretamente!")
    exit()

client = genai.Client(api_key=MINHA_CHAVE)

# --- BANCO DE DADOS DE CLIENTES (SaaS Model) ---
clientes = {
    "1": {
        "nome": "Seu Nome Teste", 
        "email": "claudinei.jb@gmail.com", 
        "especialidade": "Psiquiatria", 
        "clinica": "NeuroVida Lab", 
        "plano": "Premium", 
        "limite": 3
    },
    "2": {
        "nome": "Dr. Carlos Alberto", 
        "email": "claudinei.jb@gmail.com", 
        "especialidade": "Cardiologia", 
        "clinica": "Clinic-Ar", 
        "plano": "Basico", 
        "limite": 1
    },
    "3": {
        "nome": "Cunhada Querida", 
        "email": "claudinei.jb@gmail.com", 
        "especialidade": "Psiquiatria", 
        "clinica": "Consultorio Particular", 
        "plano": "Premium", 
        "limite": 3
    }
}

class PDF_Personalizado(FPDF):
    def __init__(self, cliente_info):
        super().__init__()
        self.cliente = cliente_info

    def header(self):
        # 1. Faixa azul do topo
        self.set_fill_color(0, 51, 102) 
        self.rect(0, 0, 210, 45, 'F') 
        
        # 2. Lógica de Logotipo Multi-Formato
        especialidade = self.cliente['especialidade']
        extensoes = [".png", ".jpg", ".jpeg"]
        logo_encontrado = None

        for ext in extensoes:
            caminhos_teste = [f"{especialidade}{ext}", f"{especialidade.lower()}{ext}"]
            for caminho in caminhos_teste:
                if os.path.exists(caminho):
                    logo_encontrado = caminho
                    break
            if logo_encontrado: break

        x_texto = 15 
        if logo_encontrado:
            self.image(logo_encontrado, x=10, y=5, w=35)
            x_texto = 55 

        # 3. Título e Subtítulo
        self.set_font("helvetica", 'B', 18)
        self.set_text_color(255, 255, 255)
        self.set_xy(x_texto, 12)
        self.cell(0, 10, "MEDICAL IN-SIGHT PREMIUM", align='L')
        
        self.ln(10)
        self.set_x(x_texto)
        self.set_font("helvetica", 'I', 10)
        self.cell(0, 5, f"Relatorio Exclusivo: {self.cliente['clinica']}", align='L')
        self.ln(25)

    def footer(self):
        self.set_y(-15)
        self.set_font("helvetica", 'I', 8)
        self.set_text_color(150, 150, 150)
        self.cell(0, 10, f"Preparado para {self.cliente['nome']} | Pagina {self.page_no()}", align='C')

def enviar_email(destinatario, nome_medico, arquivo_pdf):
    print(f"\n[E-MAIL] Preparando envio para {destinatario}...")
    try:
        msg = EmailMessage()
        msg['Subject'] = f"Sua Curadoria Científica: {nome_medico}"
        msg['From'] = EMAIL_DE
        msg['To'] = destinatario
        msg.set_content(f"Olá {nome_medico},\n\nSeu boletim científico personalizado está pronto e em anexo.\n\nAtenciosamente,\nEquipe Medical In-Sight.")

        with open(arquivo_pdf, 'rb') as f:
            msg.add_attachment(f.read(), maintype='application', subtype='pdf', filename=arquivo_pdf)

        with smtplib.SMTP_SSL('smtp.gmail.com', 465) as smtp:
            smtp.login(EMAIL_DE, SENHA_DE)
            smtp.send_message(msg)
        print("--- E-mail enviado com sucesso! ---")
    except Exception as e:
        print(f"--- Erro ao enviar e-mail: {e} ---")

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
    
    escolha = input("\nSelecione o Cliente: ")
    
    if escolha in clientes:
        user = clientes[escolha]
        artigos_brutos = buscar_pubmed(user['especialidade'], user['limite'])
        
        if artigos_brutos:
            print("[2/3] Gemini 2.0 analisando e estruturando...")
            texto_para_ia = "\n\n".join([a['texto'] for a in artigos_brutos])
            
            prompt = f"""
            Aja como um curador cientifico para {user['nome']}.
            Traduza e resuma os artigos abaixo para a especialidade {user['especialidade']}.
            Use esta estrutura para cada artigo:
            - TÍTULO: (Português)
            - RESUMO EXECUTIVO:
            - CONCLUSÃO TÉCNICA:
            Termine CADA artigo com a palavra: [FONTE]
            ARTIGOS: {texto_para_ia[:8000]}
            """
            
            response = client.models.generate_content(model="gemini-2.0-flash", contents=prompt)
            
            print("[3/3] Renderizando PDF e automatizando envios...")
            pdf = PDF_Personalizado(user)
            pdf.add_page()
            
            # Dados da entrega
            pdf.set_font("helvetica", 'B', 12)
            pdf.set_text_color(0, 51, 102)
            pdf.cell(0, 10, f"DATA: {time.strftime('%d/%m/%Y')} | DESTINATARIO: {user['nome'].upper()}", new_x="LMARGIN", new_y="NEXT")
            pdf.ln(5)
            
            partes_do_texto = response.text.split("[FONTE]")
            
            for i, parte in enumerate(partes_do_texto):
                if parte.strip():
                    pdf.set_font("helvetica", size=11)
                    pdf.set_text_color(0, 0, 0)
                    texto_seguro = parte.encode('latin-1', 'replace').decode('latin-1')
                    pdf.multi_cell(0, 8, txt=texto_seguro)
                    
                    if i < len(artigos_brutos):
                        pdf.ln(2)
                        pdf.set_font("helvetica", 'B', 10)
                        pdf.set_text_color(0, 0, 255)
                        # CENTRALIZAÇÃO CORRIGIDA: w=0 e align='C'
                        pdf.cell(0, 10, txt="--- CLIQUE AQUI PARA LER O ESTUDO COMPLETO ---", 
                                 link=artigos_brutos[i]['link'], align='C', new_x="LMARGIN", new_y="NEXT")
                        pdf.ln(10)
            
            arquivo = f"Boletim_{user['nome'].replace(' ', '_')}.pdf"
            pdf.output(arquivo)
            
            # Automação de E-mail
            enviar_email(user['email'], user['nome'], arquivo)
            
            # Automação de WhatsApp (Gera o link de redirecionamento)
            texto_wa = f"Olá {user['nome']}, seu boletim científico personalizado está pronto! Acabei de enviar para seu e-mail."
            link_wa = f"https://wa.me/?text={texto_wa.replace(' ', '%20')}"
            print(f"\n[WHATSAPP] Clique para avisar o cliente: {link_wa}")
            print(f"\n--- SUCESSO! PDF gerado: {arquivo} ---")
            
    else:
        print("Opcao invalida.")

if __name__ == "__main__":
    main()