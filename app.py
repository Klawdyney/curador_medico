import time
import os
import smtplib
import sqlite3
from email.message import EmailMessage
from dotenv import load_dotenv
from google import genai
from Bio import Entrez
from fpdf import FPDF
from fpdf.enums import XPos, YPos

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

# --- ADIÇÃO A: NOVAS FUNÇÕES DE CONTROLE DE HISTÓRICO ---

def artigo_ja_enviado(cliente_id, pubmed_id):
    """Consulta o banco para ver se o médico já recebeu este artigo."""
    conexao = sqlite3.connect('medical_insight.db')
    cursor = conexao.cursor()
    cursor.execute('SELECT 1 FROM historico_envios WHERE cliente_id = ? AND pubmed_id = ?', (cliente_id, pubmed_id))
    resultado = cursor.fetchone()
    conexao.close()
    return resultado is not None

def registrar_envio(cliente_id, pubmed_id):
    """Salva o ID do artigo no histórico do médico."""
    conexao = sqlite3.connect('medical_insight.db')
    cursor = conexao.cursor()
    cursor.execute('INSERT INTO historico_envios (cliente_id, pubmed_id) VALUES (?, ?)', (cliente_id, pubmed_id))
    conexao.commit()
    conexao.close()

# --- FUNÇÕES ORIGINAIS ---

def carregar_clientes_do_banco():
    nome_banco = 'medical_insight.db'
    try:
        conexao = sqlite3.connect(nome_banco)
        cursor = conexao.cursor()
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
        # Cabeçalho Azul Médico
        self.set_fill_color(0, 51, 102) 
        self.rect(0, 0, 210, 45, 'F') 
        
        especialidade = self.cliente['especialidade']
        esp_limpa = especialidade.split(":")[0].strip()
        
        logo_encontrado = None
        for ext in [".png", ".jpg", ".jpeg"]:
            testes = [f"{esp_limpa}{ext}", f"{esp_limpa.lower()}{ext}", f"{esp_limpa.capitalize()}{ext}"]
            for t in testes:
                if os.path.exists(t):
                    logo_encontrado = t
                    break
            if logo_encontrado: break

        x_texto = 15 
        if logo_encontrado:
            try:
                self.image(logo_encontrado, x=10, y=5, h=35)
                x_texto = 60 
            except: pass 

        self.set_font("helvetica", 'B', 18)
        self.set_text_color(255, 255, 255)
        self.set_xy(x_texto, 12)
        self.cell(0, 10, text="MEDICAL IN-SIGHT PREMIUM", align='L')
        
        self.set_font("helvetica", 'I', 11)
        self.set_xy(x_texto, 22)
        clinica_txt = self.cliente['clinica'] if self.cliente['clinica'] else esp_limpa.upper()
        self.cell(0, 5, text=f"Relatorio Exclusivo: {clinica_txt}", align='L')
        
        self.set_y(50)

    def footer(self):
        self.set_y(-15)
        self.set_font("helvetica", 'I', 8)
        self.set_text_color(150, 150, 150)
        self.cell(0, 10, text=f"Preparado para {self.cliente['nome']} | Pagina {self.page_no()}", align='C')

def buscar_pubmed(tema, limite):
    print(f"\n[1/3] Buscando {limite} evidências sobre: {tema}...")
    try:
        handle = Entrez.esearch(db="pubmed", term=tema, retmax=limite, sort="relevance")
        record = Entrez.read(handle)
        artigos = []
        for id_artigo in record.get("IdList", []):
            fetch = Entrez.efetch(db="pubmed", id=id_artigo, rettype="abstract", retmode="text")
            # --- ADIÇÃO B: INCLUINDO ID PARA CONTROLE ---
            artigos.append({
                "id": id_artigo, 
                "texto": fetch.read(), 
                "link": f"https://pubmed.ncbi.nlm.nih.gov/{id_artigo}/"
            })
        return artigos
    except Exception as e:
        print(f"[ERRO PUBMED] {e}")
        return []

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

def main():
    clientes = carregar_clientes_do_banco()
    if not clientes: return

    print("\n" + "="*40 + "\n   PORTAL MEDICAL IN-SIGHT   \n" + "="*40)
    for id_c, info in clientes.items():
        print(f"{id_c}. {info['nome']} - {info['especialidade']}")
    
    escolha = input("\nSelecione o Cliente pelo ID: ")
    if escolha in clientes:
        user = clientes[escolha]
        termo_final = f"{user['especialidade']} AND ({user['keywords']})" if user['keywords'] else user['especialidade']
        
        artigos_brutos = buscar_pubmed(termo_final, user['limite'])
        
        # --- ADIÇÃO C: LÓGICA DE FILTRAGEM DE HISTÓRICO ---
        artigos_ineditos = []
        for art in artigos_brutos:
            if not artigo_ja_enviado(escolha, art['id']):
                artigos_ineditos.append(art)
            else:
                print(f"   [INFO] Artigo {art['id']} já enviado para {user['nome']}. Pulando...")

        if artigos_ineditos:
            print(f"[2/3] Gemini analisando {len(artigos_ineditos)} novos artigos...")
            texto_para_ia = "\n\n".join([a['texto'] for a in artigos_ineditos])
            
            prompt = f"""
            Aja como um curador científico de elite para o Dr(a). {user['nome']}.
            Foco da busca: {termo_final}.

            REGRA 1 - SAUDAÇÃO (Obrigatório):
            Comece o documento com uma saudação personalizada e elegante. 
            Exemplo: "Olá, Dr(a). {user['nome']}, aqui está sua curadoria técnica sobre {user['especialidade']}..."

            REGRA 2 - ESTRUTURA DOS ARTIGOS:
            Para cada estudo, use exatamente este cabeçalho:
            **X.** [NÍVEL DE EVIDÊNCIA] - Nível (Ouro/Alto/Médio)

            Inclua: **Título:**, **Revista:**, **Resumo:**, **Conclusão Médica:** e **Nível de Evidência:**.
            Termine cada análise estritamente com a palavra [FONTE].
            Tudo em Português (Brasil).
            ARTIGOS: {texto_para_ia[:10000]}
            """
            
            response = client.models.generate_content(model="gemini-2.0-flash", contents=prompt)
            
            pdf = PDF_Personalizado(user)
            pdf.set_auto_page_break(auto=True, margin=20)
            pdf.add_page()
            
            pdf.set_font("helvetica", 'B', 11)
            pdf.set_text_color(0, 51, 102)
            tema_excl = termo_final.replace(" AND ", " + ").upper()
            pdf.cell(0, 10, text=f"DATA: {time.strftime('%d/%m/%Y')} | TEMA: {tema_excl}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            pdf.ln(5)
            
            pdf.set_text_color(0, 0, 0)
            partes = response.text.split("[FONTE]")
            for i, parte in enumerate(partes):
                if parte.strip():
                    pdf.set_font("helvetica", size=10.5)
                    texto_pdf = parte.encode('latin-1', 'replace').decode('latin-1')
                    pdf.multi_cell(0, 6, text=texto_pdf)
                    
                    if i < len(artigos_ineditos):
                        # --- REGISTRA O ARTIGO NO HISTÓRICO APÓS GERAR O PDF ---
                        registrar_envio(escolha, artigos_ineditos[i]['id'])
                        
                        pdf.set_font("helvetica", 'B', 9)
                        pdf.set_text_color(0, 102, 204)
                        pdf.cell(0, 8, text="--- CLIQUE PARA ACESSAR O ESTUDO ORIGINAL ---", 
                                 link=artigos_ineditos[i]['link'], align='C', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
                        pdf.set_text_color(0, 0, 0)
                        pdf.ln(8)
            
            arquivo = f"Boletim_{user['nome'].replace(' ', '_')}.pdf"
            pdf.output(arquivo)
            enviar_email(user['email'], user['nome'], arquivo)
            print(f"\n[SUCESSO] Curadoria inédita enviada!")
        else:
            print(f"\n[AVISO] Não há artigos novos para {user['nome']} neste tema.")
    else: print("Opção inválida.")

if __name__ == "__main__":
    main()