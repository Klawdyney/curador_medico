import sqlite3
import os
from database import get_connection
from dotenv import load_dotenv

load_dotenv()

def migrar_dados():
    print("🚀 Iniciando Migração Profissional - Fase 3...")
    
    # 1. Conexões
    try:
        local_conn = sqlite3.connect('medical_insight.db')
        local_cursor = local_conn.cursor()
        
        cloud_conn = get_connection()
        cloud_cursor = cloud_conn.cursor()
    except Exception as e:
        print(f"❌ Erro na conexão: {e}")
        return

    # 2. ATUALIZAR ESTRUTURA NA NUVEM (Schema)
    # Criamos as colunas que faltam para o sistema rodar sozinho
    print("🏗️ Ajustando estrutura da nuvem para automação...")
    colunas_novas = [
        ("status_assinatura", "TEXT DEFAULT 'ATIVO'"),
        ("dias_envio", "TEXT DEFAULT 'seg,qui'"),
        ("horario_envio", "TEXT DEFAULT '08:00'")
    ]
    
    for nome_col, tipo in colunas_novas:
        try:
            cloud_cursor.execute(f"ALTER TABLE clientes ADD COLUMN {nome_col} {tipo}")
        except:
            pass # Coluna já existe

    # 3. MIGRAR CLIENTES (Sincronizando com seus nomes: clinica e limite)
    print("📋 Sincronizando dados dos clientes...")
    local_cursor.execute("""
        SELECT nome, email, whatsapp, especialidade, clinica, keywords, plano, limite 
        FROM clientes
    """)
    clientes = local_cursor.fetchall()
    
    for c in clientes:
        query = """INSERT INTO clientes 
                   (nome, email, whatsapp, especialidade, clinica, keywords, plano, limite) 
                   VALUES (%s, %s, %s, %s, %s, %s, %s, %s) 
                   ON CONFLICT (email) DO UPDATE SET 
                   clinica = EXCLUDED.clinica, plano = EXCLUDED.plano, limite = EXCLUDED.limite"""
        cloud_cursor.execute(query, c)
    
    cloud_conn.commit()

    # 4. MIGRAR HISTÓRICO
    print("📜 Migrando histórico de envios...")
    try:
        local_cursor.execute("SELECT email_cliente, pubmed_id, titulo_artigo, link_pubmed FROM historico_envios")
        historico = local_cursor.fetchall()
        for h in historico:
            cloud_cursor.execute("""INSERT INTO historico_envios 
                                    (email_cliente, pubmed_id, titulo_artigo, link_pubmed) 
                                    VALUES (%s, %s, %s, %s) ON CONFLICT DO NOTHING""", h)
        cloud_conn.commit()
    except:
        print("🟡 Histórico já estava atualizado ou tabela não existe localmente.")

    local_conn.close()
    cloud_conn.close()
    print("\n🏁 Sistema Profissional Sincronizado! O Portal está pronto.")

if __name__ == "__main__":
    migrar_dados()