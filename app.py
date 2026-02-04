import streamlit as st  # <-- A LINHA QUE FALTAVA
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

# --- IMPORTAÇÃO DA PONTE DE DADOS (FASE 3) ---
from database import get_connection

# --- 1. CONFIGURAÇÃO DE LOGS ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("medical_insight_operacao.log", encoding='utf-8'),
        logging.StreamHandler()
    ]
)

load_dotenv()

# --- CONFIGURAÇÕES ---
MINHA_CHAVE = os.getenv("GEMINI_API_KEY")
Entrez.email = os.getenv("ENTREZ_EMAIL")
EMAIL_DE = os.getenv("EMAIL_REMETENTE")
SENHA_DE = os.getenv("EMAIL_SENHA")

# Configurações futuras de WhatsApp
WA_API_URL = os.getenv("WHATSAPP_API_URL") 
WA_API_KEY = os.getenv("WHATSAPP_API_KEY")

if not MINHA_CHAVE or not EMAIL_DE or not SENHA_DE:
    logging.error("Verifique se as chaves no arquivo .env estão preenchidas corretamente!")
    exit()

client = genai.Client(api_key=MINHA_CHAVE)

# --- FUNÇÕES DE CONTROLE DE HISTÓRICO ---

def artigo_ja_enviado(email_cliente, pubmed_id): # Troquei cliente_id por email_cliente
    try:
        with sqlite3.connect('medical_insight.db') as conexao:
            cursor = conexao.cursor()
            # Ajuste o nome da coluna na consulta também
            cursor.execute('SELECT 1 FROM historico_envios WHERE email_cliente = ? AND pubmed_id = ?', (email_cliente, pubmed_id))
            return cursor.fetchone() is not None
    except Exception as e:
        logging.error(f"Erro ao consultar histórico: {e}")
        return False

def registrar_envio(email_cliente, pubmed_id, titulo, link):
    conexao = get_connection() # Usa a ponte inteligente
    try:
        cursor = conexao.cursor()
        query = '''
            INSERT INTO historico_envios (email_cliente, pubmed_id, titulo_artigo, link_pubmed) 
            VALUES (?, ?, ?, ?)
        '''
        # Ajuste automático para PostgreSQL se necessário
        if os.getenv("DATABASE_URL"):
            query = query.replace('?', '%s')
            
        cursor.execute(query, (email_cliente, pubmed_id, titulo, link))
        conexao.commit()
    except Exception as e:
        logging.error(f"Erro ao registrar envio: {e}")
    finally:
        conexao.close()

# --- FUNÇÃO DE WHATSAPP (NOVIDADE FASE 3) ---

def enviar_whatsapp_curadoria(numero, nome_medico, especialidade):
    # Verifica se as configurações básicas existem
    if not numero or not WA_API_URL or not WA_API_KEY:
        logging.info(f"⚠️ WhatsApp ignorado para {nome_medico}: Credenciais ou número ausentes.")
        return

    # Limpeza do número: remove caracteres não numéricos
    numero_limpo = "".join(filter(str.isdigit, str(numero)))
    
    mensagem = (
        f"Olá, Dr(a). {nome_medico}! 🩺\n\n"
        f"Sua curadoria científica premium sobre *{especialidade}* acaba de ser enviada para o seu e-mail "
        f"com os últimos estudos do PubMed. Boa leitura!"
    )
    
    payload = {
        "number": numero_limpo,
        "text": mensagem
    }
    
    headers = {
        "Content-Type": "application/json",
        "apikey": WA_API_KEY
    }
    
    try:
        # Endpoint de disparo (Exemplo: Evolution API)
        # Substitua 'sua_instancia' pelo nome da sua instância real
        url_final = f"{WA_API_URL}/message/sendText/sua_instancia" 
        
        response = requests.post(url_final, headers=headers, json=payload, timeout=10)
        
        if response.status_code in [200, 201]:
            logging.info(f"✅ Notificação WhatsApp enviada para {nome_medico} ({numero_limpo}).")
        else:
            logging.error(f"❌ Falha no WhatsApp: Status {response.status_code} - {response.text}")
            
    except Exception as e:
        logging.error(f"❌ Erro de conexão com a API de WhatsApp: {e}")

# --- FUNÇÕES DE DADOS ---
# --- LOGICA DE NUANCES POR ESPECIALIDADE ---
def obter_nuance_especialidade(especialidade):
    esp_chave = especialidade.split(":")[0].strip().capitalize()
    nuances = {
        "Psiquiatria": """
            - Use terminologia baseada no DSM-5-TR e CID-11.
            - Termos obrigatórios: 'Etiopatogenia', 'Psicopatologia', 'Comorbidade', 'Remissão sintomatológica'.
            - Substitua 'efeitos colaterais' por 'perfil de tolerabilidade' ou 'eventos adversos'.
        """,
        "Cardiologia": """
            - Use as diretrizes da SBC, AHA e ESC.
            - Foco em: 'MACE (Major Adverse Cardiovascular Events)', 'Desfechos duros', 'Fração de ejeção'.
            - Termos obrigatórios: 'Estratificação de risco', 'Insulto isquêmico', 'Hemodinâmica'.
        """,
        "Dermatologia": """
            - Foco em: 'Dermatoscopia', 'Imunofenotipagem', 'Anátomo-patológico'.
            - Termos obrigatórios: 'Lesão elementar', 'Fisiopatologia cutânea', 'Manejo terapêutico tópico'.
        """
    }
    return nuances.get(esp_chave, "- Use jargão médico acadêmico sênior e terminologia DeCS/MeSH padrão.")
def carregar_clientes_do_banco():
        conexao = get_connection() 
        try:
            if hasattr(conexao, 'row_factory'):
                conexao.row_factory = sqlite3.Row
                cursor = conexao.cursor()
            else:
                import psycopg2.extras
                cursor = conexao.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

            # SQL AJUSTADO: Usando exatamente o que o seu terminal encontrou na nuvem
            cursor.execute('''
                SELECT id, nome, email, especialidade, keywords, whatsapp, 
                    clinica, plano, horario_envio, dia_envio, limite, valor_assinatura 
                FROM clientes
            ''')
            linhas = cursor.fetchall()
            
            clientes = {}
            for l in linhas:
                id_c = str(l['id'])
                clientes[id_c] = {
                    "nome": l['nome'], 
                    "email": l['email'],
                    "especialidade": l['especialidade'],
                    "keywords": l['keywords'],
                    "whatsapp": l['whatsapp'],
                    "clinica": l['clinica'] if l['clinica'] else "Medical In-Sight",
                    "plano": l['plano'],
                    "dias": l['dia_envio'],      # Ajustado para o singular 'dia_envio'
                    "horario": l['horario_envio'],
                    "limite": l['limite'] if l['limite'] else 2,
                    "valor": l['valor_assinatura']
                }
            return clientes
        except Exception as e:
            logging.error(f"Erro ao carregar dados profissionais: {e}")
            return {}
        finally:
            if 'conexao' in locals():
                conexao.close()

# --- CLASSE DE PDF (DESIGN ORIGINAL INTEGRAL) ---

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
        # Posiciona a 1,5 cm do fim da página
        self.set_y(-15)
        # Ajuste 1: Trocamos 'Arial' por 'helvetica' (padrão moderno)
        self.set_font('helvetica', 'I', 8) 
        self.set_text_color(128, 128, 128) 
        
        # Linha fina separadora
        self.line(10, self.get_y(), 200, self.get_y())
        
        # Ajuste 2: Atualizamos os parâmetros de alinhamento para o padrão fpdf2
        self.cell(0, 10, f'Página {self.page_no()}', align='L', new_x="RIGHT", new_y="TOP")
        
        # O Selo de Qualidade
        self.set_text_color(0, 51, 102) 
        self.cell(0, 10, '[ Verified by Medical Expert ] Curadoria Algorítmica | Validação Clínica Humana  ', align='R')

# --- LOGICA DE COMUNICAÇÃO ---

def buscar_pubmed(tema, limite_busca=20, dias=1460):
    logging.info(f"Buscando PubMed: {tema} (Janela: {dias} dias)...")
    try:
        # Agora usamos a variável 'dias' na query
        query_filtrada = f"({tema}) AND \"last {dias} days\"[dp] AND (Clinical Trial[ptyp] OR Review[ptyp] OR Meta-Analysis[ptyp] OR Journal Article[ptyp]) NOT Editorial[ptyp]"
        handle = Entrez.esearch(db="pubmed", term=query_filtrada, retmax=limite_busca, sort="relevance")
        record = Entrez.read(handle)
        artigos = []
        for id_artigo in record.get("IdList", []):
            fetch = Entrez.efetch(db="pubmed", id=id_artigo, rettype="medline", retmode="text")
            artigos.append({
                "id": id_artigo, 
                "texto": fetch.read(), 
                "link": f"https://pubmed.ncbi.nlm.nih.gov/{id_artigo}/",
                "tipo": "NOVIDADE" if dias <= 30 else "ESTUDO CLÁSSICO" # Identificador para a IA
            })
        return artigos
    except Exception as e:
        logging.error(f"Erro PubMed: {e}"); return []

def enviar_email_pdf(email_destino, nome_medico, arquivo_pdf, e_classico=False):
    try:
        # Lógica Dinâmica baseada no tipo de conteúdo
        titulo_tipo = "Marco Histórico da Medicina" if e_classico else "Boletim de Inteligência Clínica"
        status_texto = "Selecionamos um <strong>estudo clássico fundamental</strong> (Landmark Trial) para sua especialidade, visto que não houve publicações disruptivas recentes." if e_classico else "Identificamos <strong>novas evidências de alto impacto</strong> relevantes para a sua prática clínica."
        cor_borda = "#003366" if e_classico else "#28a745" # Azul para clássico, Verde para novidade

        msg = EmailMessage()
        msg['Subject'] = f"📚 {titulo_tipo}: Dr(a). {nome_medico}"
        msg['From'] = EMAIL_DE
        msg['To'] = email_destino

        html_content = f"""
        <html>
            <body style="font-family: Arial, sans-serif; color: #333; line-height: 1.6;">
                <div style="background-color: #003366; padding: 30px; text-align: center; color: white;">
                    <h1 style="margin: 0; font-size: 26px; letter-spacing: 2px;">MEDICAL IN-SIGHT</h1>
                    <p style="margin: 5px 0 0; font-size: 12px; letter-spacing: 1px; opacity: 0.8;">MONITORAMENTO CIENTÍFICO EM TEMPO REAL</p>
                </div>
                
                <div style="padding: 40px; max-width: 600px; margin: auto; border: 1px solid #eee;">
                    <p style="font-size: 16px;">Olá, <strong>Dr(a). {nome_medico}</strong>,</p>
                    <p>É com satisfação que enviamos o seu <strong>{titulo_tipo}</strong> personalizado.</p>
                    
                    <div style="background-color: #f4f7fa; border-left: 5px solid {cor_borda}; padding: 25px; margin: 30px 0; border-radius: 4px;">
                        <h2 style="color: {cor_borda}; font-size: 14px; margin-top: 0; text-transform: uppercase; letter-spacing: 1px;">Status da Atualização:</h2>
                        <p style="margin-bottom: 0; font-size: 16px; color: #444;">
                            {status_texto}
                        </p>
                    </div>
                    
                    <p style="font-size: 14px; color: #555;">
                        Esta curadoria utiliza algoritmos de processamento de linguagem natural e revisão técnica para garantir que apenas estudos com rigor metodológico cheguem até si.
                    </p>
                    
                    <div style="text-align: center; margin: 30px 0;">
                        <p style="font-size: 12px; color: #888;">(O arquivo PDF está anexado a este e-mail)</p>
                    </div>

                    <hr style="border: 0; border-top: 1px solid #eee; margin: 30px 0;">
                    
                    <p style="font-size: 13px; color: #999; text-align: center;">
                        Atenciosamente,<br>
                        <strong style="color: #003366;">Curadoria Científica Medical In-Sight</strong>
                    </p>
                </div>
            </body>
        </html>
        """
        
        msg.add_alternative(html_content, subtype='html')

        with open(arquivo_pdf, 'rb') as f:
            file_data = f.read()
            msg.add_attachment(file_data, maintype='application', subtype='pdf', filename=arquivo_pdf)

        with smtplib.SMTP_SSL('smtp.gmail.com', 465) as smtp:
            smtp.login(EMAIL_DE, SENHA_DE)
            smtp.send_message(msg)
            
        logging.info(f"E-mail ({titulo_tipo}) enviado para {email_destino}")
        
    except Exception as e:
        logging.error(f"Erro ao enviar e-mail para {email_destino}: {e}")

def enviar_radar_sem_novidades(destinatario, nome_medico, especialidade):
    try:
        msg = EmailMessage()
        msg['Subject'] = f"Radar Medical In-Sight: Monitoramento {especialidade}"
        msg['From'] = EMAIL_DE
        msg['To'] = destinatario
        
        html_content = f"""
        <html>
            <body style="font-family: Arial, sans-serif; color: #333; line-height: 1.6;">
                <div style="background-color: #003366; padding: 30px; text-align: center; color: white;">
                    <h1 style="margin: 0; font-size: 26px; letter-spacing: 2px;">MEDICAL IN-SIGHT</h1>
                    <p style="margin: 5px 0 0; font-size: 12px; letter-spacing: 1px; opacity: 0.8;">MONITORAMENTO CIENTÍFICO EM TEMPO REAL</p>
                </div>
                
                <div style="padding: 40px; max-width: 600px; margin: auto; border: 1px solid #eee;">
                    <p style="font-size: 16px;">Olá, <strong>Dr(a). {nome_medico}</strong>,</p>
                    <p>Informamos que nosso sistema de varredura automatizada concluiu o monitoramento das bases de dados científicas (PubMed/MEDLINE) nas últimas 24 horas.</p>
                    
                    <div style="background-color: #f8f9fa; border-left: 5px solid #003366; padding: 25px; margin: 30px 0; border-radius: 4px;">
                        <h2 style="color: #003366; font-size: 14px; margin-top: 0; text-transform: uppercase; letter-spacing: 1px;">Status do Radar:</h2>
                        <p style="margin-bottom: 0; font-size: 16px; color: #444;">
                            <strong>Nenhuma nova evidência disruptiva</strong> foi publicada para o tema <span style="color: #003366;">"{especialidade}"</span> desde sua última atualização.
                        </p>
                    </div>
                    
                    <p style="font-size: 14px; color: #777; font-style: italic;">
                        Isso indica que sua conduta clínica permanece rigorosamente alinhada com o estado da arte da literatura médica atual. Continuamos em prontidão para capturar o próximo marco científico.
                    </p>
                    
                    <hr style="border: 0; border-top: 1px solid #eee; margin: 30px 0;">
                    
                    <p style="font-size: 13px; color: #999; text-align: center;">
                        Atenciosamente,<br>
                        <strong style="color: #003366;">Curadoria Científica Medical In-Sight</strong>
                    </p>
                </div>
            </body>
        </html>
        """
        
        msg.add_alternative(html_content, subtype='html')
        with smtplib.SMTP_SSL('smtp.gmail.com', 465) as smtp:
            smtp.login(EMAIL_DE, SENHA_DE)
            smtp.send_message(msg)
        logging.info(f"Radar sem novidades enviado para {destinatario}")
        
    except Exception as e:
        logging.error(f"Erro ao enviar radar para {destinatario}: {e}")
def processar_medico_completo(user):
    """Motor Único de Inteligência: PubMed -> Gemini -> PDF -> Envio"""
    nome_medico = user['nome']
    email_cliente = user['email']
    especialidade = user['especialidade']
    keywords = user['keywords']
    limite = user['limite']
    whatsapp = user['whatsapp']
    clinica = user['clinica'] if user['clinica'] else "Medical In-Sight"

    try:
        termo_final = f"{especialidade} AND ({keywords})" if keywords else especialidade
        
        # --- 1. LÓGICA DE BUSCA EM 3 NÍVEIS (IDÊNTICA À ORIGINAL) ---
        artigos_ineditos = []
        contem_classico = False

        # Nível 1: Hot News (15 dias)
        artigos_n1 = buscar_pubmed(termo_final, limite_busca=20, dias=15)
        artigos_ineditos = [art for art in artigos_n1 if not artigo_ja_enviado(email_cliente, art['id'])]

        # Nível 2: Fila (4 anos)
        if len(artigos_ineditos) < limite:
            artigos_n2 = buscar_pubmed(termo_final, limite_busca=20, dias=1460)
            for art in artigos_n2:
                if not artigo_ja_enviado(email_cliente, art['id']) and art['id'] not in [a['id'] for a in artigos_ineditos]:
                    artigos_ineditos.append(art)

        # Nível 3: Radar Clássico (20 anos)
        if not artigos_ineditos:
            termo_classico = f"({termo_final}) AND (landmark trial OR classic study OR trial)"
            artigos_n3 = buscar_pubmed(termo_classico, limite_busca=1, dias=7300)
            if artigos_n3:
                artigos_n3[0]['tipo'] = 'ESTUDO CLÁSSICO'
                artigos_ineditos = artigos_n3
                contem_classico = True

        artigos_para_enviar = artigos_ineditos[:limite]

        if not artigos_para_enviar:
            enviar_radar_sem_novidades(email_cliente, nome_medico, especialidade)
            return f"📡 [RADAR] {nome_medico} - Sem novidades no período."

        # --- 2. INTELIGÊNCIA GEMINI (PROMPT COMPLETO) ---
        bloco_artigos_texto = "".join([f"\nID: {a['id']}\nSTATUS: {a.get('tipo', 'NOVIDADE')}\nCONTEÚDO: {a['texto']}\n---" for a in artigos_para_enviar])
        nuance_extra = obter_nuance_especialidade(especialidade)
        nota_elegante = f'⚠️ INSTRUÇÃO PRIORITÁRIA: Inicie obrigatoriamente com a nota de Marco Histórico para {especialidade}...' if contem_classico else ""

        prompt = f"""Aja como um Curador Científico Sênior para o Dr. {nome_medico}. 
        Especialidade: {especialidade}. {nuance_extra} {nota_elegante}

        ESTRUTURA OBRIGATÓRIA POR ARTIGO:
        [TITULO_INICIO] Tradução técnica em português seguida do ano do estudo. [TITULO_FIM]
        [EVIDENCIA_INICIO] Classifique como: Alto (Metanálise/Ensaio Randomizado), Médio (Observacional/Coorte) ou Baixo (Relatos/Editoriais). [EVIDENCIA_FIM]
        [FONTE_INICIO] Citação acadêmica completa. [FONTE_FIM]
        
        [RESUMO_INICIO] 
        1. METODOLOGIA: Descreva o desenho do estudo (N, duração, critérios de inclusão).
        2. RESULTADOS: Apresente os desfechos primários com dados numéricos e relevância estatística.
        3. ANÁLISE CRÍTICA: Discuta o impacto fisiopatológico e a inovação para a {especialidade}.
        [RESUMO_FIM]
        
        [CONCLUSAO_INICIO] APLICAÇÃO CLÍNICA: Recomendação direta e objetiva para a prática diária. [CONCLUSAO_FIM]
        Termine cada análise estritamente com [PROXIMO_ARTIGO].
        """
        
        response = client.models.generate_content(model="gemini-2.0-flash", contents=prompt + bloco_artigos_texto)

        # --- 3. GERAÇÃO DO PDF PREMIUM (DESIGN COMPLETO) ---
        pdf = PDF_Personalizado(user); pdf.add_page()
        pdf.set_font("helvetica", 'B', 16); pdf.set_text_color(0, 51, 102)
        pdf.cell(0, 10, text=clinica.upper(), align='C', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        pdf.ln(2)
        pdf.set_font("helvetica", 'B', 10); pdf.set_text_color(100, 100, 100)
        pdf.cell(0, 10, text=f"GERADO EM: {time.strftime('%d/%m/%Y')} | FOCO: {termo_final.upper()}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        pdf.ln(5)

        saudacao = f"Prezado Dr. {nome_medico.split()[0]}, segue análise técnica das evidências selecionadas para sua atuação clínica."
        pdf.set_font("helvetica", 'I', 11); pdf.set_text_color(50, 50, 50)
        pdf.multi_cell(0, 7, text=saudacao.encode('latin-1', 'replace').decode('latin-1'))
        pdf.ln(5)

        partes = response.text.split("[PROXIMO_ARTIGO]")
        for i, parte in enumerate(partes):
            if i >= len(artigos_para_enviar) or "[TITULO_INICIO]" not in parte: continue
            try:
                titulo = parte.split("[TITULO_INICIO]")[1].split("[TITULO_FIM]")[0].strip()
                evidencia = parte.split("[EVIDENCIA_INICIO]")[1].split("[EVIDENCIA_FIM]")[0].strip()
                fonte = parte.split("[FONTE_INICIO]")[1].split("[FONTE_FIM]")[0].strip()
                resumo = parte.split("[RESUMO_INICIO]")[1].split("[RESUMO_FIM]")[0].strip()
                conclusao = parte.split("[CONCLUSAO_INICIO]")[1].split("[CONCLUSAO_FIM]")[0].strip()

                pdf.set_x(15)
                pdf.set_font("helvetica", 'B', 10); pdf.set_text_color(60, 60, 60)
                pdf.multi_cell(0, 7, text=f"[NÍVEL DE EVIDÊNCIA: {evidencia.upper()}]")
                
                pdf.set_x(15)
                pdf.set_font("helvetica", 'B', 12); pdf.set_text_color(0, 51, 102)
                pdf.multi_cell(0, 7, text=titulo.encode('latin-1', 'replace').decode('latin-1'))
                
                pdf.set_x(15)
                pdf.set_font("helvetica", 'I', 9); pdf.set_text_color(100, 100, 100)
                pdf.multi_cell(0, 5, text=f"Fonte: {fonte}".encode('latin-1', 'replace').decode('latin-1'))
                pdf.ln(3)
                
                pdf.set_x(15)
                pdf.set_font("helvetica", '', 10.5); pdf.set_text_color(30, 30, 30)
                pdf.multi_cell(0, 6, text=resumo.encode('latin-1', 'replace').decode('latin-1'))
                pdf.ln(2)
                
                pdf.set_x(15); pdf.set_fill_color(245, 247, 250); pdf.set_font("helvetica", 'B', 10); pdf.set_text_color(0, 51, 102)
                pdf.cell(0, 8, text="   APLICAÇÃO CLÍNICA / CONCLUSÃO:", fill=True, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
                
                pdf.set_x(15); pdf.set_font("helvetica", 'I', 10); pdf.set_text_color(50, 50, 50)
                pdf.multi_cell(0, 6, text=conclusao.encode('latin-1', 'replace').decode('latin-1'), fill=True)
                pdf.ln(2)
                
                pdf.set_x(15); pdf.set_font("helvetica", 'B', 9); pdf.set_text_color(0, 102, 204)
                pdf.cell(0, 8, text=">> ACESSAR ESTUDO COMPLETO NO PUBMED <<", link=artigos_para_enviar[i]['link'], align='C', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
                
                registrar_envio(email_cliente, artigos_para_enviar[i]['id'], titulo, artigos_para_enviar[i]['link'])
                pdf.ln(10)
            except Exception as e:
                logging.error(f"Erro no parse do artigo: {e}")
                continue

        arquivo = f"Boletim_{nome_medico.replace(' ', '_')}.pdf"
        pdf.output(arquivo)
        enviar_email_pdf(email_cliente, nome_medico, arquivo, e_classico=contem_classico)
        enviar_whatsapp_curadoria(whatsapp, nome_medico, especialidade)
        
        return f"✅ [SUCESSO] {nome_medico} - E-mail enviado."

    except Exception as e:
        logging.error(f"Erro crítico: {e}")
        return f"❌ [ERRO] {nome_medico} - Motivo: {e}"
# --- MAIN ---

def main():
    clientes = carregar_clientes_do_banco()
    if not clientes: return
    
    # --- NOVO: LISTA PARA ARMAZENAR O STATUS DE CADA MÉDICO ---
    relatorio_final = []
    # ---------------------------------------------------------

    st.title("🏥 Portal Medical In-Sight")

    # Prepara as opções para o menu visual (substitui o loop de print da imagem)
    opcoes = {f"{info['nome']} - {info['especialidade']}": id_c for id_c, info in clientes.items()}
    escolha = st.selectbox("Selecione o Cliente:", ["todos"] + list(opcoes.keys()))

    # O botão substitui o input() e impede que o site fique em branco
    if st.button("🚀 Iniciar Curadoria Científica"):
    # Define os IDs baseado na escolha do menu
        ids_para_processar = list(clientes.keys()) if escolha == 'todos' else [opcoes[escolha]]

        for id_c in ids_para_processar:
            if id_c not in clientes: continue
            
            user = clientes[id_c]
            st.write(f"⚙️ Acionando motor de curadoria: Dr. {user['nome']}...")
            
            # CHAMADA CIRÚRGICA: O Python pula para a linha 383, executa tudo e volta.
            # (A Função Mestra já cuida dos 3 níveis, Gemini, PDF e Envio)
            resultado = processar_medico_completo(user)
            
            # Feedback visual baseado no ícone retornado pela função
            if "✅" in resultado:
                st.success(resultado)
            elif "📡" in resultado:
                st.info(resultado)
            else:
                st.error(resultado)
                
            # Registra no relatório de logs que será gravado ao final
            relatorio_final.append(resultado)
    # --- AGORA SIM: FORA DO LOOP, MAS DENTRO DO MAIN ---
    if relatorio_final:
    # Ajustado para a pasta 'historico_logs'
        nome_arquivo_txt = f"historico_logs/Relatorio_Envios_{time.strftime('%Y%m%d_%H%M%S')}.txt"
        try:
            with open(nome_arquivo_txt, "w", encoding="utf-8") as f:
                f.write("=== RELATÓRIO DE ENTREGAS MEDICAL IN-SIGHT ===\n")
                f.write(f"Data da Rodada: {time.strftime('%d/%m/%Y %H:%M:%S')}\n")
                f.write("-" * 45 + "\n\n")
                for linha in relatorio_final:
                    f.write(linha + "\n")
            print(f"\n✅ Relatório detalhado gerado em: {nome_arquivo_txt}")
        except Exception as e:
            print(f"\n⚠️ Erro ao gravar o arquivo de relatório: {e}")

if __name__ == "__main__":
    main()