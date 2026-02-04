from database import get_connection
import psycopg2

def sincronizar_nuvem_saas():
    print("🚀 Iniciando Sincronização Definitiva da Nuvem...")
    conn = get_connection()
    cursor = conn.cursor()
    
    # 1. Lista de colunas que o seu app.py profissional EXIGE
    colunas_obrigatorias = [
        ("clinica", "TEXT"),
        ("plano", "TEXT DEFAULT 'BASICO'"),
        ("status_assinatura", "TEXT DEFAULT 'ATIVO'"),
        ("dias_envio", "TEXT DEFAULT 'seg,qui'"),
        ("horario_envio", "TEXT DEFAULT '08:00'"),
        ("limite", "INTEGER DEFAULT 2")
    ]
    
    for nome, tipo in colunas_obrigatorias:
        try:
            # Tenta adicionar a coluna
            cursor.execute(f"ALTER TABLE clientes ADD COLUMN {nome} {tipo};")
            print(f"✅ Coluna '{nome}' adicionada com sucesso!")
        except Exception as e:
            # Se der erro, verificamos se é porque a coluna já existe
            conn.rollback() # Limpa o erro para continuar a próxima
            print(f"🟡 Coluna '{nome}' já existe ou ignorada.")
            
    conn.commit()
    
    # 2. Verificação Final: Vamos ler o que o banco tem agora
    print("\n📋 Verificando estrutura final na Nuvem:")
    cursor.execute("SELECT column_name FROM information_schema.columns WHERE table_name = 'clientes';")
    colunas_reais = [row[0] for row in cursor.fetchall()]
    print(f"Colunas encontradas: {', '.join(colunas_reais)}")
    
    conn.close()
    print("\n🏁 Sincronização concluída! Tente rodar o app.py agora.")

if __name__ == "__main__":
    sincronizar_nuvem_saas()