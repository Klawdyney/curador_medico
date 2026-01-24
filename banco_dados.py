import sqlite3

def criar_banco():
    conexao = sqlite3.connect('medical_insight.db')
    cursor = conexao.cursor()
    
    # Criando a tabela de clientes
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS clientes (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        nome TEXT NOT NULL,
        email TEXT NOT NULL,
        especialidade TEXT NOT NULL,
        clinica TEXT,
        plano TEXT DEFAULT 'Gratuito',
        limite INTEGER DEFAULT 1
    )
    ''')
    
    conexao.commit()
    conexao.close()
    print("Banco de dados criado com sucesso!")

def adicionar_cliente(nome, email, especialidade, clinica, plano, limite):
    conexao = sqlite3.connect('medical_insight.db')
    cursor = conexao.cursor()
    cursor.execute('''
    INSERT INTO clientes (nome, email, especialidade, clinica, plano, limite)
    VALUES (?, ?, ?, ?, ?, ?)
    ''', (nome, email, especialidade, clinica, plano, limite))
    
    conexao.commit()
    conexao.close()
    print(f"Cliente {nome} cadastrado com sucesso!")

if __name__ == "__main__":
    criar_banco()
    # Exemplo de cadastro inicial (pode apagar ou alterar depois)
    adicionar_cliente("Cunhada Querida", "claudinei.jb@gmail.com", "Psiquiatria", "Consultorio Particular", "Premium", 5)
    adicionar_cliente("Estudante Teste", "claudinei.jb@gmail.com", "Cardiologia", "UFPR", "Estudante", 2)