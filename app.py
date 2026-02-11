import time
import os
import smtplib
import sqlite3
import re
import logging
import requests 
import resend # Adicione este import no topo do seu app.py
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
resend.api_key = os.getenv("RESEND_API_KEY")

# Configurações futuras de WhatsApp
WA_API_URL = os.getenv("WHATSAPP_API_URL") 
WA_API_KEY = os.getenv("WHATSAPP_API_KEY")
# Validação Profissional: Agora checa o Resend em vez do Gmail
if not MINHA_CHAVE or not resend.api_key:
    logging.error("ERRO: Verifique se GEMINI_API_KEY e RESEND_API_KEY estão no seu .env!")
    exit()

client = genai.Client(api_key=MINHA_CHAVE)

# --- FUNÇÕES DE CONTROLE DE HISTÓRICO ---

# --- FUNÇÕES DE CONTROLE DE HISTÓRICO (VERSÃO NUVEM SUPABASE) ---

def artigo_ja_enviado(email_cliente, pubmed_id):
    """Verifica no SUPABASE se o artigo já foi enviado."""
    conexao = get_connection() # Agora usa a conexão da nuvem, não o arquivo local
    try:
        cursor = conexao.cursor()
        
        # Query ajustada para PostgreSQL (Supabase usa %s)
        query = "SELECT 1 FROM historico_envios WHERE email_cliente = %s AND pubmed_id = %s"
        cursor.execute(query, (email_cliente, str(pubmed_id)))
        
        return cursor.fetchone() is not None
    except Exception as e:
        logging.error(f"⚠️ Erro ao consultar histórico no Supabase: {e}")
        return False # Na dúvida, envia (melhor pecar pelo excesso do que pela falta)
    finally:
        conexao.close()

def registrar_envio(email_cliente, pubmed_id, titulo, link):
    """Grava no SUPABASE que o envio foi feito."""
    conexao = get_connection()
    try:
        cursor = conexao.cursor()
        
        # Timestamp atual para auditoria
        data_hoje = time.strftime('%Y-%m-%d %H:%M:%S')
        
        # Inserção correta no PostgreSQL
        query = '''
            INSERT INTO historico_envios (email_cliente, pubmed_id, titulo_artigo, link_pubmed, data_envio) 
            VALUES (%s, %s, %s, %s, %s)
        '''
        cursor.execute(query, (email_cliente, str(pubmed_id), titulo, link, data_hoje))
        conexao.commit()
        logging.info(f"💾 Histórico salvo no Supabase: {pubmed_id}")
    except Exception as e:
        logging.error(f"❌ Erro ao registrar envio no Supabase: {e}")
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
        """,
        "Oftalmologia": """
            - Foco em: 'OCT (Tomografia de Coerência Óptica)', 'Pressão Intraocular (PIO)', 'Aguidade Visual'.
            - Termos obrigatórios: 'Fundo de olho', 'Biomicroscopia', 'Segmento anterior/posterior'.
            - Priorize estudos sobre degeneração macular, glaucoma e retinopatias.
        """,
        "Cirurgia geral": """
            - Foco em: 'Técnicas minimamente invasivas', 'Morbimortalidade pós-operatória', 'Videolaparoscopia'.
            - Termos obrigatórios: 'Manejo cirúrgico', 'Complicações pós-operatórias', 'Tempo cirúrgico'.
            - Destaque avanços em cirurgia robótica e trauma.
        """,
        "Ginecologia e obstetrícia": """
            - Foco em: 'Desfechos perinatais', 'Vitalidade fetal', 'Assistência pré-natal'.
            - Termos obrigatórios: 'Ciclo gravídico-puerperal', 'Hormonioterapia', 'Medicina Fetal'.
            - Priorize diretrizes da FEBRASGO e ACOG.
        """,
        "Radiologia": """
            - Foco em: 'Sensibilidade e Especificidade', 'Laudo estruturado', 'Achados incidentais'.
            - Termos obrigatórios: 'Critérios BI-RADS/LI-RADS/TI-RADS', 'Meio de contraste', 'Modalidade de imagem (RM/TC/US)'.
            - Enfatize a relevância clínica de novos protocolos de imagem.
        """,
        "Anestesiologia": """
            - Foco em: 'Manejo de via aérea', 'Estabilidade hemodinâmica perioperatória', 'Analgesia multimodal'.
            - Termos obrigatórios: 'Recuperação pós-anestésica (RPA)', 'Bloqueio neuroaxial', 'Farmacocinética dos anestésicos'.
            - Priorize segurança do paciente e protocolos ERAS.
        """,
        "Oncologia": """
            - Use critérios RECIST para resposta tumoral.
            - Foco em: 'Sobrevida Global (OS)', 'Sobrevida Livre de Progressão (PFS)', 'Imunoterápicos'.
            - Termos obrigatórios: 'Estadiamento TNM', 'Mutação driver', 'Terapia-alvo'.
        """,
        "Neurologia": """
            - Foco em: 'Neuroimagem', 'Déficit focal', 'Escala NIHSS'.
            - Termos obrigatórios: 'Neuroplasticidade', 'Fisiopatologia sináptica', 'Padrão eletroencefalográfico'.
        """,
        "Pediatria": """
            - Foco em: 'Marcos do desenvolvimento', 'Faixa etária pediátrica'.
            - Use doses baseadas em mg/kg quando aplicável e referências de Puericultura.
        """,
        "Ortopedia": """
            - Foco em: 'Biomecânica', 'Consolidação óssea', 'Cinética'.
            - Termos obrigatórios: 'Manejo cirúrgico', 'Redução anatômica', 'Osteossíntese'.
        """,
        "Padrao": """
            - Use rigor científico absoluto e terminologia médica acadêmica de alto nível.
            - Extraia obrigatoriamente: Metodologia, Tamanho da Amostra (N), P-valor e Intervalo de Confiança.
            - Mantenha o foco estrito em Medicina Baseada em Evidências (MBE).
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
        # --- 1. BLINDAGEM JURÍDICA (Fica acima da sua linha) ---
        # Posiciona 2.5 cm do fim para caber o aviso legal
        self.set_y(-25) 
        self.set_font('helvetica', '', 6) # Letra bem pequena e discreta
        self.set_text_color(160, 160, 160) # Cinza bem claro para não poluir
        disclaimer = "ISENÇÃO DE RESPONSABILIDADE: Conteúdo gerado por IA para fins estritamente informativos. Não substitui diretrizes oficiais ou julgamento clínico. Verifique as fontes originais."
        # O multi_cell quebra o texto se for longo e centraliza
        self.multi_cell(0, 3, text=disclaimer.encode('latin-1', 'replace').decode('latin-1'), align='C')
        
        # --- 2. O SEU DESIGN ORIGINAL (Mantido Intacto) ---
        # Posiciona na altura padrão que você gosta (-1.5 cm)
        self.set_y(-15)
    
        # Linha fina separadora
        self.set_draw_color(200, 200, 200) # Garante a cor cinza da linha
        self.line(10, self.get_y(), 200, self.get_y())
        
        # Configuração da fonte do rodapé
        self.set_font('helvetica', 'I', 8) 
        self.set_text_color(128, 128, 128) 
        
        # Página (Seu código original)
        self.cell(0, 10, f'Página {self.page_no()}', align='L', new_x="RIGHT", new_y="TOP")
        
        # O Selo de Qualidade (Seu texto completo original)
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
        # 1. LÓGICA DINÂMICA ORIGINAL PRESERVADA
        titulo_tipo = "Marco Histórico da Medicina" if e_classico else "Boletim de Inteligência Clínica"
        status_texto = "Selecionamos um <strong>estudo clássico fundamental</strong> (Landmark Trial) para sua especialidade, visto que não houve publicações disruptivas recentes." if e_classico else "Identificamos <strong>novas evidências de alto impacto</strong> relevantes para a sua prática clínica."
        cor_borda = "#003366" if e_classico else "#28a745"

        # 2. SEU HTML INTEGRAL (Design Premium Mantido)
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
                    <div style="margin-top: 30px; padding-top: 20px; border-top: 1px solid #eee; font-size: 10px; color: #aaa; text-align: center;">
                        <p>DISCLAIMER: Este e-mail contém resumos gerados por Inteligência Artificial. A decisão clínica é soberana do médico.</p>
                        <p><a href="#" style="color: #aaa;">Descadastrar (Unsubscribe)</a> - <a href="#" style="color: #aaa;">Termos de Uso</a></p>
                    </div>
                </div>
            </body>
        </html>
        """

        # 3. PREPARAÇÃO DO ANEXO (Formato exigido pelo Resend)
        with open(arquivo_pdf, 'rb') as f:
            file_data = f.read()

        # 4. DISPARO VIA API (Substitui o smtplib.SMTP_SSL que causava erro)
        resend.Emails.send({
            "from": "Medical In-Sight <curadoria@medinsight.com.br>",
            "to": [email_destino],
            "subject": f"📚 {titulo_tipo}: Dr(a). {nome_medico}",
            "html": html_content,
            "attachments": [
                {
                    "filename": arquivo_pdf,
                    "content": list(file_data) # O Resend Python SDK exige lista de bytes
                }
            ]
        })
            
        logging.info(f"✅ Sucesso Premium: E-mail ({titulo_tipo}) enviado via Resend para {email_destino}")
        
    except Exception as e:
        logging.error(f"❌ Erro ao enviar e-mail via Resend para {email_destino}: {e}")

def enviar_radar_sem_novidades(destinatario, nome_medico, especialidade):
    try:
        # SEU HTML INTEGRAL E VALIDADO (Design Premium Mantido)
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
                    <div style="margin-top: 30px; padding-top: 20px; border-top: 1px solid #eee; font-size: 10px; color: #aaa; text-align: center;">
                        <p>DISCLAIMER: Este e-mail contém resumos gerados por Inteligência Artificial. A decisão clínica é soberana do médico.</p>
                        <p><a href="#" style="color: #aaa;">Descadastrar (Unsubscribe)</a> - <a href="#" style="color: #aaa;">Termos de Uso</a></p>
                    </div>
                </div>
            </body>
        </html>
        """
        
        # DISPARO VIA API RESEND (Substituindo smtplib e login do Gmail)
        resend.Emails.send({
            "from": "Medical In-Sight <curadoria@medinsight.com.br>",
            "to": [destinatario],
            "subject": f"Radar Medical In-Sight: Monitoramento {especialidade}",
            "html": html_content
        })
        
        logging.info(f"✅ Radar sem novidades enviado via Resend para {destinatario}")
        
    except Exception as e:
        logging.error(f"❌ Erro ao enviar radar via Resend para {destinatario}: {e}")

def traduzir_para_ingles_medico(termo_pt):
    """Usa a IA para converter termos (PT -> EN) com blindagem dupla de modelos e retries."""
    if not termo_pt: return ""
    
    # Lista de modelos para fugir de bloqueios regionais ou erros 404
    modelos_para_tentar = ["gemini-1.5-flash", "gemini-2.0-flash"]
    
    for modelo_atual in modelos_para_tentar:
        for tentativa in range(2): # Tenta 2 vezes cada modelo
            try:
                prompt = f"Translate this medical term from Portuguese to English (MeSH term) for PubMed search. Output ONLY the English term: {termo_pt}"
                
                # Chamada com o modelo da vez
                response = client.models.generate_content(model=modelo_atual, contents=prompt)
                termo_en = response.text.strip()
                
                logging.info(f"🌐 Tradução Concluída ({modelo_atual}): '{termo_pt}' -> '{termo_en}'")
                return termo_en
                
            except Exception as e:
                # Se for erro de localização ou modelo não encontrado, ele avisa e tentará o próximo
                logging.warning(f"⚠️ Falha no modelo {modelo_atual} (Tentativa {tentativa+1}): {e}")
                time.sleep(2) 
                
    # Se todos os modelos falharem, retorna o original para a busca não parar
    logging.error(f"❌ Todos os modelos de tradução falharam para: {termo_pt}")
    return termo_pt
    
def processar_medico_completo(user):
    """Motor Único de Inteligência: PubMed -> Gemini -> PDF -> Envio"""
    nome_medico = user['nome']
    email_cliente = user['email']
    especialidade = user['especialidade']
    keywords_pt = user['keywords'] # Guardamos o original em PT
    
    # 1. Correção do Limite (Segurança para não dar erro se vier texto)
    try:
        limite = int(user['limite'])
    except:
        limite = 2

    whatsapp = user['whatsapp']
    clinica = user['clinica'] if user['clinica'] else "Medical In-Sight"

    try:
        # 2. Tradução Invisível (PT -> EN)
        # O robô traduz o que o médico digitou para o inglês antes de buscar
        keywords_en = traduzir_para_ingles_medico(keywords_pt)
        especialidade_en = traduzir_para_ingles_medico(especialidade)

        # O termo_final agora fica em INGLÊS para o PubMed encontrar tudo
        termo_final = f"{especialidade_en} AND ({keywords_en})" if keywords_en else especialidade_en
        
        logging.info(f"🔎 Busca Oficial (EN): {termo_final} | Médico digitou: {keywords_pt}")
        
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
        
        # --- 2. INTELIGÊNCIA GEMINI COM BLINDAGEM DUPLA (FLASH 1.5 e 2.0) ---
        response = None
        # Lista de modelos para tentar fugir de bloqueios regionais ou 404
        modelos_para_tentar = ["gemini-1.5-flash", "gemini-2.0-flash"]
        
        for modelo_atual in modelos_para_tentar:
            if response: break # Se já conseguiu, sai do loop de modelos
            
            for tentativa in range(2): # Tenta 2 vezes cada modelo
                try:
                    logging.info(f"🤖 Tentando análise com o modelo: {modelo_atual}")
                    response = client.models.generate_content(
                        model=modelo_atual, 
                        contents=prompt + bloco_artigos_texto
                    )
                    break # Sucesso! Sai do loop de tentativas
                except Exception as e:
                    logging.warning(f"⚠️ Erro no modelo {modelo_atual}: {e}")
                    time.sleep(3) # Pausa curta antes de tentar o próximo

        if not response:
            return f"❌ [ERRO] {nome_medico} - Todos os modelos de IA falharam (Localização ou Limite)."

        # --- 3. GERAÇÃO DO PDF PREMIUM (SEGUE O CÓDIGO ABAIXO) ---

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

# def main():
#     clientes = carregar_clientes_do_banco()
#     if not clientes: return
#     
#     st.title("🏥 Portal Medical In-Sight")
#     # ... (as linhas comentadas continuam aqui)

# if __name__ == "__main__":
#     main()