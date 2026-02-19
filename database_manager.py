import psycopg2
import psycopg2.extras # Importante para transformar dados do banco em 'dicionários'
import os
from dotenv import load_dotenv

load_dotenv()

def get_connection():
    """Estabelece a conexão segura com o Supabase na Nuvem."""
    url = os.getenv("DATABASE_URL")
    if not url:
        raise ValueError("❌ Erro: DATABASE_URL não encontrada no .env!")
    return psycopg2.connect(url, sslmode='require')

# --- FUNÇÕES PARA O PORTAL (CLIENTE) ---

def cadastrar_medico(dados):
    """Realiza o primeiro cadastro do médico no sistema."""
    query = """
    INSERT INTO clientes (nome, email, whatsapp, especialidade, clinica, keywords, plano, limite, dia_envio, horario_envio, senha)
    VALUES (%(nome)s, %(email)s, %(whatsapp)s, %(especialidade)s, %(clinica)s, %(keywords)s, %(plano)s, %(limite)s, %(dia_envio)s, %(horario_envio)s, %(senha)s)
    """
    try:
        with get_connection() as conn:
            with conn.cursor() as cur:
                if 'senha' not in dados: dados['senha'] = '123456'
                cur.execute(query, dados)
                conn.commit()
                return True
    except Exception as e:
        print(f"Erro ao cadastrar: {e}")
        return False

def buscar_medico_por_email(email):
    """Busca os dados de um médico específico para carregar no Painel de Perfil."""
    query = "SELECT * FROM clientes WHERE email = %s"
    try:
        with get_connection() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute(query, (email,))
                return cur.fetchone()
    except Exception as e:
        print(f"Erro ao buscar perfil: {e}")
        return None

def atualizar_perfil_medico(email_original, novos_dados):
    """
    Salva TODAS as alterações que o médico fez no seu perfil.
    """
    query = """
    UPDATE clientes 
    SET plano = %(plano)s, 
        limite = %(limite)s, 
        especialidade = %(especialidade)s, -- Adicionada a peça que faltava!
        keywords = %(keywords)s, 
        dia_envio = %(dia_envio)s, 
        horario_envio = %(horario_envio)s,
        clinica = %(clinica)s
    WHERE email = %(email_original)s
    """
    try:
        novos_dados['email_original'] = email_original
        with get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(query, novos_dados)
                conn.commit()
                return True
    except Exception as e:
        print(f"Erro técnico ao atualizar perfil: {e}")
        return False

# --- FUNÇÕES PARA O ROBÔ (WORKER) ---

def buscar_todos_os_medicos_ativos():
    """Retorna a lista de médicos para o processamento do robô worker.py."""
    query = "SELECT * FROM clientes WHERE ativo = TRUE"
    try:
        with get_connection() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute(query)
                return cur.fetchall()
    except Exception as e:
        print(f"Erro ao buscar médicos ativos: {e}")
        return []