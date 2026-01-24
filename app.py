import time
import os
import smtplib
import sqlite3
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

# --- FUNÇÃO PARA LER O BANCO DE DADOS ---
def carregar_clientes_do_banco():
    nome_banco = 'medical_insight.db'
    try:
        conexao = sqlite3.connect(nome_banco)
        cursor = conexao.cursor()
        # Captura todas as colunas necessárias, incluindo a nova 'keywords'
        cursor.execute('SELECT id, nome, email, especialidade, clinica, plano, limite, keywords FROM clientes')
        linhas = cursor.fetchall()
        conexao.close()
        
        return {str(l[0]): {
            "nome": l[1], 
            "email": l[2], 
            "especialidade": l[3], 
            "clinica": l[4], 
            "plano": l[5], 
            "limite": l[6],
            "keywords": l[7] 
        } for l in linhas}
    except Exception as e:
        print(f"\n[ERRO] Falha ao ler banco de dados: {e}")
        return {}

class PDF_Personalizado(FPDF):
    def __init__(self, cliente_info):
        super().__init__()
        self.cliente = cliente_info

    def header(self):
        self.set_fill_color(0, 51, 102) 
        self.rect(0, 0, 210, 45, 'F') 
        
        especialidade = self.cliente['especialidade']
        especialidade_limpa = especialidade.split(":")[0].strip()
        
        extensoes = [".png", ".jpg", ".jpeg"]
        logo_encontrado = None

        for ext in extensoes:
            caminhos_teste = [f"{especialidade_limpa}{ext}", f"{especialidade_limpa.lower()}{ext}"]
            for caminho in caminhos_teste:
                if os.path.exists(caminho):
                    logo_encontrado = caminho
                    break
            if logo_encontrado: break

        x_texto = 15 
        if logo_encontrado:
            self.image(logo_encontrado, x=10, y=5, w=35)
            x_texto = 55 

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
    print(f"\n[E-MAIL] Enviando para {destinatario}...")
    try:
        msg = EmailMessage()
        msg['Subject'] = f"Sua Curadoria Científica: {nome_medico}"
        msg['From'] = EMAIL_DE
        msg['To'] = destinatario
        msg.set_content(f"Olá {nome_medico},\n\nSeu boletim personalizado está em anexo.")
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
    try:
        handle = Entrez.esearch(db="pubmed", term=tema, retmax=limite, sort="relevance")
        record = Entrez.read(handle)
        artigos = []
        for id_artigo in record.get("IdList", []):
            fetch = Entrez.efetch(db="pubmed", id=id_artigo, rettype="abstract", retmode="text")
            artigos.append({"texto": fetch.read(), "link": f"https://pubmed.ncbi.nlm.nih.gov/{id_artigo}/"})
        return artigos
    except Exception as e:
        print(f"[ERRO PUBMED] Falha na busca: {e}")
        return []

def main():
    clientes = carregar_clientes_do_banco()
    if not clientes:
        print("\n[ERRO] Nenhum cliente encontrado no banco de dados!")
        return

    print("\n" + "="*40 + "\n   PORTAL MEDICAL IN-SIGHT (DB MODE)   \n" + "="*40)
    for id_c, info in clientes.items():
        print(f"{id_c}. {info['nome']} - {info['especialidade']}")
    
    escolha = input("\nSelecione o Cliente pelo ID: ")
    if escolha in clientes:
        user = clientes[escolha]
        
        # TRATAMENTO DE VALOR NULO (NULL) VINDO DO BANCO
        adicional = user.get('keywords')
        if adicional is None:
            adicional = ""
        
        termo_final = user['especialidade']
        
        # Se houver palavras-chave, combina com a especialidade
        if adicional.strip():
            termo_final = f"{user['especialidade']} AND ({adicional})"
        
        artigos_brutos = buscar_pubmed(termo_final, user['limite'])
        
        if artigos_brutos:
            print("[2/3] Gemini 2.0 analisando e traduzindo...")
            texto_para_ia = "\n\n".join([a['texto'] for a in artigos_brutos])
            
            prompt = f"""
            Aja como curador científico de elite para {user['nome']}. 
            Foco da busca: {termo_final}.

            INSTRUÇÕES OBRIGATÓRIAS:
            1. TRADUZA TUDO para o Português (Brasil), inclusive os TÍTULOS.
            2. Se houver termos específicos ({adicional}), priorize evidências sobre eles no resumo.
            3. Estrutura por artigo: TÍTULO, RESUMO, CONCLUSÃO. 
            4. Termine cada análise com [FONTE].

            ARTIGOS: {texto_para_ia[:8000]}
            """
            
            response = client.models.generate_content(model="gemini-2.0-flash", contents=prompt)
            
            print("[3/3] Gerando PDF e enviando...")
            pdf = PDF_Personalizado(user)
            pdf.add_page()
            
            # Formatação do Título do PDF (Removendo caracteres técnicos de busca)
            tema_exibicao = termo_final.replace(" AND ", " E ").replace("(", "").replace(")", "").upper()
            pdf.set_font("helvetica", 'B', 12)
            pdf.cell(0, 10, f"DATA: {time.strftime('%d/%m/%Y')} | TEMA: {tema_exibicao}", ln=1)
            pdf.ln(5)
            
            partes = response.text.split("[FONTE]")
            for i, parte in enumerate(partes):
                if parte.strip():
                    pdf.set_font("helvetica", size=11)
                    # Codificação segura para caracteres latinos
                    texto_pdf = parte.encode('latin-1', 'replace').decode('latin-1')
                    pdf.multi_cell(0, 8, txt=texto_pdf)
                    
                    if i < len(artigos_brutos):
                        pdf.set_font("helvetica", 'B', 10)
                        pdf.set_text_color(0, 0, 255)
                        pdf.cell(0, 10, txt="--- CLIQUE AQUI PARA LER O ESTUDO ORIGINAL (INGLES) ---", 
                                 link=artigos_brutos[i]['link'], align='C', ln=1)
                        pdf.set_text_color(0, 0, 0)
                        pdf.ln(10)
            
            arquivo = f"Boletim_{user['nome'].replace(' ', '_')}.pdf"
            pdf.output(arquivo)
            enviar_email(user['email'], user['nome'], arquivo)
            
            print(f"\n[SUCESSO] Curadoria sobre '{tema_exibicao}' enviada!")
        else:
            print(f"\n[AVISO] Nenhum artigo encontrado para: {termo_final}")
    else:
        print("Opção inválida.")

if __name__ == "__main__":
    main()