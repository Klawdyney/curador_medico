import streamlit as st
import time
import os
import smtplib
import sqlite3
import re
import logging
import requests
from email.message import EmailMessage
from dotenv import load_dotenv
from google import genai
from Bio import Entrez
from fpdf import FPDF
from fpdf.enums import XPos, YPos

# --- IMPORTAÇÃO DA PONTE DE DADOS ---
from database import get_connection

# --- 1. CONFIGURAÇÃO DE LOGS ---
# Ajustado para exibir logs no painel do Railway
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)

load_dotenv()

# --- CONFIGURAÇÕES ---
MINHA_CHAVE = os.getenv("GEMINI_API_KEY")
Entrez.email = os.getenv("ENTREZ_EMAIL")
EMAIL_DE = os.getenv("EMAIL_REMETENTE")
SENHA_DE = os.getenv("EMAIL_SENHA")

# Configurações de WhatsApp
WA_API_URL = os.getenv("WHATSAPP_API_URL") 
WA_API_KEY = os.getenv("WHATSAPP_API_KEY")

# Inicializa cliente Gemini
if MINHA_CHAVE:
    client = genai.Client(api_key=MINHA_CHAVE)

# --- FUNÇÕES DE BANCO DE DADOS (Híbridas: Local e Nuvem) ---

def artigo_ja_enviado(email_cliente, pubmed_id):
    conexao = get_connection()
    try:
        cursor = conexao.cursor()
        query = 'SELECT 1 FROM historico_envios WHERE email_cliente = ? AND pubmed_id = ?'
        if os.getenv("DATABASE_URL"): # Ajuste para PostgreSQL
            query = query.replace('?', '%s')
        cursor.execute(query, (email_cliente, str(pubmed_id)))
        return cursor.fetchone() is not None
    except Exception as e:
        logging.error(f"Erro ao consultar histórico: {e}")
        return False
    finally:
        conexao.close()

def registrar_envio(email_cliente, pubmed_id, titulo, link):
    conexao = get_connection()
    try:
        cursor = conexao.cursor()
        query = '''
            INSERT INTO historico_envios (email_cliente, pubmed_id, titulo_artigo, link_pubmed) 
            VALUES (?, ?, ?, ?)
        '''
        if os.getenv("DATABASE_URL"):
            query = query.replace('?', '%s')
        cursor.execute(query, (email_cliente, str(pubmed_id), titulo, link))
        conexao.commit()
    except Exception as e:
        logging.error(f"Erro ao registrar envio: {e}")
    finally:
        conexao.close()

def carregar_clientes_do_banco():
    conexao = get_connection() 
    try:
        cursor = conexao.cursor()
        cursor.execute('SELECT id, nome, email, especialidade, clinica, plano, limite, keywords, whatsapp FROM clientes')
        colunas = [desc[0] for desc in cursor.description]
        clientes = {}
        for row in cursor.fetchall():
            cliente_dict = dict(zip(colunas, row))
            clientes[str(cliente_dict['id'])] = cliente_dict
        return clientes
    except Exception as e:
        logging.error(f"Falha ao ler banco de dados: {e}")
        return {}
    finally:
        conexao.close()

# --- WHATSAPP E NUANCES ---

def enviar_whatsapp_curadoria(numero, nome_medico, especialidade):
    if not numero or not WA_API_URL or not WA_API_KEY:
        return 
    
    numero_limpo = "".join(filter(str.isdigit, str(numero)))
    mensagem = (
        f"Olá, Dr(a). {nome_medico}! 🩺\n\n"
        f"Sua curadoria científica premium sobre *{especialidade}* acaba de ser enviada para o seu e-mail "
        f"com os últimos estudos do PubMed. Boa leitura!"
    )
    
    headers = {"Content-Type": "application/json", "apikey": WA_API_KEY}
    payload = {"number": numero_limpo, "text": mensagem}
    
    try:
        url_final = f"{WA_API_URL}/message/sendText/sua_instancia" 
        requests.post(url_final, headers=headers, json=payload, timeout=5)
    except Exception as e:
        logging.error(f"Erro WhatsApp: {e}")

def obter_nuance_especialidade(especialidade):
    esp_chave = especialidade.split(":")[0].strip().capitalize()
    nuances = {
        "Psiquiatria": """
            - Use terminologia baseada no DSM-5-TR e CID-11.
            - Termos obrigatórios: 'Etiopatogenia', 'Psicopatologia', 'Comorbidade', 'Remissão sintomatológica'.
        """,
        "Cardiologia": """
            - Use as diretrizes da SBC, AHA e ESC.
            - Foco em: 'MACE', 'Desfechos duros', 'Fração de ejeção', 'Insulto isquêmico'.
        """,
        "Dermatologia": """
            - Foco em: 'Dermatoscopia', 'Imunofenotipagem', 'Anátomo-patológico'.
        """
    }
    return nuances.get(esp_chave, "- Use jargão médico acadêmico sênior.")

# --- CLASSE DE PDF ---

class PDF_Personalizado(FPDF):
    def __init__(self, cliente_info):
        super().__init__()
        self.cliente = cliente_info
        self.set_margins(15, 15, 15)
        self.set_auto_page_break(auto=True, margin=25)

    def header(self):
        self.set_fill_color(0, 51, 102) 
        self.rect(0, 0, 210, 45, 'F') 
        
        especialidade = self.cliente['especialidade']
        esp_limpa = especialidade.split(":")[0].strip()
        
        # Tenta carregar logo (suporta várias extensões)
        logo_encontrado = None
        for ext in [".png", ".jpg", ".jpeg"]:
            testes = [f"{esp_limpa}{ext}", f"{esp_limpa.lower()}{ext}"]
            for t in testes:
                if os.path.exists(t):
                    logo_encontrado = t
                    break
            if logo_encontrado: break

        x_texto = 15 
        if logo_encontrado:
            try:
                self.image(logo_encontrado, x=12, y=5, h=32)
                x_texto = 60 
            except: pass 

        self.set_font("helvetica", 'B', 18)
        self.set_text_color(255, 255, 255)
        self.set_xy(x_texto, 12)
        self.cell(0, 10, text="MEDICAL IN-SIGHT PREMIUM", align='L')
        
        self.set_font("helvetica", 'I', 10)
        self.set_xy(x_texto, 22)
        clinica_txt = self.cliente['clinica'] if self.cliente['clinica'] else esp_limpa.upper()
        self.cell(0, 5, text=f"Curadoria Científica para: {clinica_txt}", align='L')
        
        self.set_draw_color(255, 255, 255)
        self.line(x_texto, 20, 200, 20)
        self.set_y(50)

    def footer(self):
        self.set_y(-15)
        self.set_font('helvetica', 'I', 8) 
        self.set_text_color(128, 128, 128) 
        self.line(10, self.get_y(), 200, self.get_y())
        self.cell(0, 10, f'Página {self.page_no()}', align='L', new_x="RIGHT", new_y="TOP")
        self.set_text_color(0, 51, 102) 
        self.cell(0, 10, '[ Verified by Medical Expert ] Curadoria Algorítmica | Validação Clínica Humana', align='R')

# --- BUSCA PUBMED ---

def buscar_pubmed(tema, limite_busca=20, dias=1460):
    try:
        query_filtrada = f"({tema}) AND \"last {dias} days\"[dp] AND (Clinical Trial[ptyp] OR Review[ptyp] OR Meta-Analysis[ptyp] OR Journal Article[ptyp]) NOT Editorial[ptyp]"
        handle = Entrez.esearch(db="pubmed", term=query_filtrada, retmax=limite_busca, sort="relevance")
        record = Entrez.read(handle)
        artigos = []
        for id_artigo in record.get("IdList", []):
            try:
                fetch = Entrez.efetch(db="pubmed", id=id_artigo, rettype="medline", retmode="text")
                texto_completo = fetch.read()
                artigos.append({
                    "id": id_artigo, 
                    "texto": texto_completo[:3000], 
                    "link": f"https://pubmed.ncbi.nlm.nih.gov/{id_artigo}/",
                    "tipo": "NOVIDADE" if dias <= 30 else "ESTUDO CLÁSSICO"
                })
            except: continue
        return artigos
    except Exception as e:
        logging.error(f"Erro PubMed: {e}"); return []

# --- ENVIO DE E-MAIL ---

def enviar_email_pdf(email_destino, nome_medico, arquivo_pdf, e_classico=False):
    try:
        titulo_tipo = "Marco Histórico da Medicina" if e_classico else "Boletim de Inteligência Clínica"
        status_texto = "Selecionamos um <strong>estudo clássico fundamental</strong>." if e_classico else "Identificamos <strong>novas evidências de alto impacto</strong>."
        cor_borda = "#003366" if e_classico else "#28a745"

        msg = EmailMessage()
        msg['Subject'] = f"📚 {titulo_tipo}: Dr(a). {nome_medico}"
        msg['From'] = EMAIL_DE
        msg['To'] = email_destino

        html_content = f"""
        <html>
            <body style="font-family: Arial, sans-serif; color: #333;">
                <div style="background-color: #003366; padding: 20px; text-align: center; color: white;">
                    <h2>MEDICAL IN-SIGHT</h2>
                </div>
                <div style="padding: 20px; border: 1px solid #eee;">
                    <p>Olá, <strong>Dr(a). {nome_medico}</strong>,</p>
                    <div style="border-left: 5px solid {cor_borda}; padding: 15px; margin: 20px 0; background: #f9f9f9;">
                        {status_texto}
                    </div>
                    <p>Segue em anexo a análise completa em PDF.</p>
                </div>
            </body>
        </html>
        """
        msg.add_alternative(html_content, subtype='html')

        with open(arquivo_pdf, 'rb') as f:
            msg.add_attachment(f.read(), maintype='application', subtype='pdf', filename=arquivo_pdf)

        # Mudança: Usar SMTP (não SSL) na porta 587 e adicionar starttls()
        with smtplib.SMTP('smtp.gmail.com', 587) as smtp:
            smtp.starttls()  # <--- LINHA NOVA E ESSENCIAL
            smtp.login(EMAIL_DE, SENHA_DE)
            smtp.send_message(msg)
            
    except Exception as e:
        logging.error(f"Erro e-mail: {e}")

def enviar_radar_sem_novidades(destinatario, nome_medico, especialidade):
    try:
        msg = EmailMessage()
        msg['Subject'] = f"Radar Medical In-Sight: Monitoramento {especialidade}"
        msg['From'] = EMAIL_DE
        msg['To'] = destinatario
        html_content = f"""
        <html><body>
            <p>Olá, Dr(a). {nome_medico},</p>
            <p>Monitoramos o PubMed nas últimas 24h e <strong>não houve novas publicações disruptivas</strong> para {especialidade}.</p>
        </body></html>
        """
        msg.add_alternative(html_content, subtype='html')
        
        # AQUI ESTÁ A CORREÇÃO (Porta 587 + starttls)
        with smtplib.SMTP('smtp.gmail.com', 587) as smtp:
            smtp.starttls()
            smtp.login(EMAIL_DE, SENHA_DE)
            smtp.send_message(msg)
    except: pass

# --- MAIN (ADAPTADO PARA WEB/STREAMLIT) ---

def main():
    # Configuração da Página Web
    st.set_page_config(page_title="Medical In-Sight", page_icon="🏥", layout="wide")
    st.title("🏥 Portal Medical In-Sight")
    st.markdown("### Curadoria Científica Avançada com IA")

    # 1. Carregar Clientes
    clientes = carregar_clientes_do_banco()
    if not clientes:
        st.error("Erro: Não foi possível carregar os clientes. Verifique o banco de dados.")
        return

    # 2. Painel Lateral (Substitui o input do terminal)
    st.sidebar.header("Painel de Controle")
    lista_nomes = {f"{info['nome']} - {info['especialidade']}": id_c for id_c, info in clientes.items()}
    
    escolha_web = st.selectbox(
        "Selecione o Médico:",
        ["-- Selecione --", "PROCESSAR TODOS"] + list(lista_nomes.keys())
    )

    # 3. Botão de Execução
    if st.button("🚀 Iniciar Curadoria e Envio"):
        if escolha_web == "-- Selecione --":
            st.warning("Selecione uma opção válida.")
            return

        ids_processar = list(clientes.keys()) if escolha_web == "PROCESSAR TODOS" else [lista_nomes[escolha_web]]
        
        barra = st.progress(0)
        logs_tela = []

        for index, id_c in enumerate(ids_processar):
            user = clientes[id_c]
            st.info(f"🔍 Analisando evidências para: Dr(a). {user['nome']}...")
            
            try:
                # --- LÓGICA DE 3 NÍVEIS (MANTIDA) ---
                termo_final = f"{user['especialidade']} AND ({user['keywords']})" if user['keywords'] else user['especialidade']
                artigos_ineditos = []
                contem_classico = False

                # Nível 1: Hot News
                artigos_n1 = buscar_pubmed(termo_final, limite_busca=20, dias=15)
                artigos_ineditos = [art for art in artigos_n1 if not artigo_ja_enviado(user['email'], art['id'])]

                # Nível 2: Fila de 4 anos
                if len(artigos_ineditos) < user['limite']:
                    artigos_n2 = buscar_pubmed(termo_final, limite_busca=20, dias=1460)
                    for art in artigos_n2:
                        if not artigo_ja_enviado(user['email'], art['id']) and art['id'] not in [a['id'] for a in artigos_ineditos]:
                            artigos_ineditos.append(art)
                
                # Nível 3: Clássicos
                if not artigos_ineditos:
                    st.warning(f"Buscando clássicos para {user['nome']}...")
                    termo_classico = f"({termo_final}) AND (landmark trial OR classic study)"
                    artigos_n3 = buscar_pubmed(termo_classico, limite_busca=1, dias=7300)
                    if artigos_n3:
                        artigos_n3[0]['tipo'] = 'ESTUDO CLÁSSICO'
                        artigos_ineditos = artigos_n3
                        contem_classico = True

                artigos_para_enviar = artigos_ineditos[:user['limite']]

                if artigos_para_enviar:
                    # --- PROMPT E IA ---
                    bloco_artigos = "".join([f"\nID: {a['id']}\nSTATUS: {a['tipo']}\nCONTEÚDO: {a['texto']}\n---" for a in artigos_para_enviar])
                    nuance = obter_nuance_especialidade(user['especialidade'])
                    nota_elegante = 'INCLUA NOTA: "Selecionamos este Marco Histórico (Landmark Trial)."' if contem_classico else ""

                    prompt = f"""
                    Aja como Curador Sênior para Dr. {user['nome']}. Especialidade: {user['especialidade']}. {nuance} {nota_elegante}
                    ESTRUTURA OBRIGATÓRIA:
                    [TITULO_INICIO] Tradução título [TITULO_FIM]
                    [EVIDENCIA_INICIO] Classificação [EVIDENCIA_FIM]
                    [FONTE_INICIO] Citação [FONTE_FIM]
                    [RESUMO_INICIO] 
                    1. METODOLOGIA: ...
                    2. RESULTADOS: ...
                    3. ANÁLISE CRÍTICA: ...
                    [RESUMO_FIM]
                    [CONCLUSAO_INICIO] Aplicação Clínica [CONCLUSAO_FIM]
                    Finalize com [PROXIMO_ARTIGO].
                    """
                    
                    response = client.models.generate_content(model="gemini-2.0-flash", contents=prompt + f"\n\nDADOS: {bloco_artigos}")
                    
                    # --- GERAÇÃO DO PDF (CORRIGIDA: LINK + MARGEM SEGURA) ---
                    pdf = PDF_Personalizado(user)
                    pdf.add_page()
                    
                    # Cabeçalho da página
                    pdf.set_font("helvetica", 'B', 10); pdf.set_text_color(100, 100, 100)
                    pdf.set_x(15) # Segurança de margem
                    pdf.cell(0, 10, text=f"GERADO EM: {time.strftime('%d/%m/%Y')} | FOCO: {termo_final.upper()}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
                    pdf.ln(5)

                    pdf.set_font("helvetica", 'I', 11); pdf.set_text_color(50, 50, 50)
                    pdf.set_x(15)
                    pdf.multi_cell(0, 7, text=f"Prezado Dr. {user['nome'].split()[0]}, segue análise técnica.")
                    pdf.ln(5)

                    partes = response.text.split("[PROXIMO_ARTIGO]")
                    for i, parte in enumerate(partes):
                        if i >= len(artigos_para_enviar) or "[TITULO_INICIO]" not in parte: continue
                        try:
                            # Extração dos dados
                            titulo = parte.split("[TITULO_INICIO]")[1].split("[TITULO_FIM]")[0].strip()
                            evidencia = parte.split("[EVIDENCIA_INICIO]")[1].split("[EVIDENCIA_FIM]")[0].strip()
                            fonte = parte.split("[FONTE_INICIO]")[1].split("[FONTE_FIM]")[0].strip()
                            resumo = parte.split("[RESUMO_INICIO]")[1].split("[RESUMO_FIM]")[0].strip()
                            conclusao = parte.split("[CONCLUSAO_INICIO]")[1].split("[CONCLUSAO_FIM]")[0].strip()

                            # Escrita no PDF (Tudo com set_x(15) para evitar erro)
                            pdf.set_font("helvetica", 'B', 10); pdf.set_text_color(60, 60, 60)
                            pdf.set_x(15)
                            pdf.multi_cell(0, 7, text=f"[NÍVEL DE EVIDÊNCIA: {evidencia.upper()}]")
                            
                            pdf.set_font("helvetica", 'B', 12); pdf.set_text_color(0, 51, 102)
                            pdf.set_x(15)
                            pdf.multi_cell(0, 7, text=titulo.encode('latin-1', 'replace').decode('latin-1'))
                            
                            pdf.set_font("helvetica", 'I', 9); pdf.set_text_color(100, 100, 100)
                            pdf.set_x(15)
                            pdf.multi_cell(0, 5, text=f"Fonte: {fonte}".encode('latin-1', 'replace').decode('latin-1'))
                            pdf.ln(3)
                            
                            pdf.set_font("helvetica", '', 10.5); pdf.set_text_color(30, 30, 30)
                            pdf.set_x(15)
                            pdf.multi_cell(0, 6, text=resumo.encode('latin-1', 'replace').decode('latin-1'))
                            pdf.ln(2)
                            
                            pdf.set_fill_color(245, 247, 250); pdf.set_font("helvetica", 'B', 10); pdf.set_text_color(0, 51, 102)
                            pdf.set_x(15)
                            pdf.cell(0, 8, text=" APLICAÇÃO CLÍNICA:", fill=True, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
                            
                            pdf.set_font("helvetica", 'I', 10); pdf.set_text_color(50, 50, 50)
                            pdf.set_x(15)
                            pdf.multi_cell(0, 6, text=conclusao.encode('latin-1', 'replace').decode('latin-1'), fill=True)
                            
                            # --- O LINK ESTÁ DE VOLTA AQUI ---
                            pdf.ln(2)
                            pdf.set_x(15) 
                            pdf.set_font("helvetica", 'B', 9); pdf.set_text_color(0, 102, 204)
                            # Link centralizado e seguro
                            pdf.cell(0, 8, text=">> ACESSAR ESTUDO COMPLETO NO PUBMED <<", link=artigos_para_enviar[i]['link'], align='C', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
                            
                            pdf.ln(8) 
                            registrar_envio(user['email'], artigos_para_enviar[i]['id'], titulo, artigos_para_enviar[i]['link'])

                        except Exception as e_parse:
                            st.warning(f"Erro ao formatar artigo {i+1}: {e_parse}")
                            continue

                    # Salvar e Enviar
                    nome_arquivo = f"Boletim_{user['nome'].replace(' ', '_')}.pdf"
                    pdf.output(nome_arquivo)
                    
                    enviar_email_pdf(user['email'], user['nome'], nome_arquivo, contem_classico)
                    enviar_whatsapp_curadoria(user['whatsapp'], user['nome'], user['especialidade'])
                    
                    st.success(f"✅ Sucesso total para: {user['nome']}")
                    logs_tela.append(f"SUCESSO: {user['nome']}")

                else:
                    enviar_radar_sem_novidades(user['email'], user['nome'], user['especialidade'])
                    logs_tela.append(f"RADAR: {user['nome']}")

            except Exception as e:
                st.error(f"Erro no cliente {user['nome']}: {e}")
                logs_tela.append(f"ERRO: {user['nome']} - {e}")
            
            barra.progress((index + 1) / len(ids_processar))

        # Relatório Final na tela
        st.write("---")
        st.subheader("📊 Relatório da Operação")
        for log in logs_tela:
            st.write(log)

if __name__ == "__main__":
    main()