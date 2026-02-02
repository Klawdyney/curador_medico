import sqlite3
import os
from database import get_connection
from dotenv import load_dotenv

load_dotenv()

def migrar_dados():
    print("🚀 Iniciando migração para a Fase 3...")
    
    # 1. Conecta ao Banco Local (SQLite)
    local_conn = sqlite3.connect('medical_insight.db')
    local_cursor = local_conn.cursor()
    
    # 2. Conecta ao Banco Nuvem (PostgreSQL)
    try:
        cloud_conn = get_connection()
        if not os.getenv("DATABASE_URL"):
            print("❌ ERRO: DATABASE_URL não encontrada no .env.")
            return
        cloud_cursor = cloud_conn.cursor()
    except Exception as e:
        print(f"❌ Erro ao conectar na nuvem: {e}")
        return

    # --- PASSO 3: MIGRAR CLIENTES (Seleção Explícita) ---
    print("📋 Migrando clientes...")
    # Selecionamos exatamente as 8 colunas que vamos inserir
    local_cursor.execute("""
        SELECT nome, email, whatsapp, especialidade, clinica, keywords, plano, limite 
        FROM clientes
    """)
    clientes = local_cursor.fetchall()
    
    for c in clientes:
        # Usamos exatamente 8 marcadores %s para as 8 colunas
        query = """INSERT INTO clientes 
                   (nome, email, whatsapp, especialidade, clinica, keywords, plano, limite) 
                   VALUES (%s, %s, %s, %s, %s, %s, %s, %s) 
                   ON CONFLICT (email) DO NOTHING"""
        cloud_cursor.execute(query, c)
    
    cloud_conn.commit()
    print(f"✅ {len(clientes)} clientes processados.")

    # --- PASSO 4: MIGRAR HISTÓRICO (Essencial para não repetir artigos) ---
    print("📜 Migrando histórico de envios...")
    local_cursor.execute("""
        SELECT email_cliente, pubmed_id, titulo_artigo, link_pubmed 
        FROM historico_envios
    """)
    historico = local_cursor.fetchall()
    
    for h in historico:
        query_h = """INSERT INTO historico_envios 
                     (email_cliente, pubmed_id, titulo_artigo, link_pubmed) 
                     VALUES (%s, %s, %s, %s)"""
        cloud_cursor.execute(query_h, h)
    
    cloud_conn.commit()
    print(f"✅ {len(historico)} registros de histórico migrados.")

    # --- FINALIZAÇÃO ---
    local_conn.close()
    cloud_conn.close()
    print("\n🏁 Migração concluída com sucesso! Seu produto está 100% na nuvem.")

if __name__ == "__main__":
    migrar_dados()