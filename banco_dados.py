import sqlite3

def criar_banco():
    """
    Cria a estrutura do banco de dados Medical In-Sight.
    Garante que todas as colunas necessárias, incluindo 'keywords', existam.
    """
    conexao = sqlite3.connect('medical_insight.db')
    cursor = conexao.cursor()
    
    # Criando a tabela de clientes com a estrutura completa
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS clientes (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        nome TEXT NOT NULL,
        email TEXT NOT NULL,
        especialidade TEXT NOT NULL,
        clinica TEXT,
        plano TEXT DEFAULT 'Básico',
        limite INTEGER DEFAULT 1,
        keywords TEXT
    )
    ''')
    
    conexao.commit()
    conexao.close()
    print("✅ Estrutura do banco de dados verificada/criada com sucesso!")

def adicionar_cliente_teste(nome, email, especialidade, clinica, plano, limite, keywords):
    """
    Função utilitária para inserir um perfil de teste rapidamente.
    """
    try:
        conexao = sqlite3.connect('medical_insight.db')
        cursor = conexao.cursor()
        cursor.execute('''
            INSERT INTO clientes (nome, email, especialidade, clinica, plano, limite, keywords)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', (nome, email, especialidade, clinica, plano, limite, keywords))
        conexao.commit()
        conexao.close()
        print(f"🚀 Perfil de teste '{nome}' adicionado com sucesso!")
    except Exception as e:
        print(f"❌ Erro ao adicionar teste: {e}")

if __name__ == "__main__":
    # 1. Garante que a "casa" está pronta
    criar_banco()
    
    # 2. Cadastro de Demonstração para o Gabriel (Residente)
    # Sugestão de keywords para um residente de Clínica Médica/Neurologia
    adicionar_cliente_teste(
        nome="Gabriel Residente", 
        email="seu_email_aqui@gmail.com", # Altere para seu e-mail de teste
        especialidade="Neurologia", 
        clinica="Hospital de Clínicas", 
        plano="Premium", 
        limite=3, 
        keywords="Stroke, Multiple Sclerosis, Neuroplasticity"
    )
    
    print("\nPronto! O banco está configurado e com um perfil de teste para sua apresentação.")