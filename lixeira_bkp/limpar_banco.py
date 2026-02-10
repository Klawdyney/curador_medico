import sqlite3

def limpeza_total():
    conexao = sqlite3.connect('medical_insight.db')
    cursor = conexao.cursor()
    
    # Apaga todos os clientes e todo o histórico de e-mails já enviados
    cursor.execute("DELETE FROM clientes")
    cursor.execute("DELETE FROM historico_envios")
    
    # Reseta os IDs para o Gabriel ser o ID 1
    cursor.execute("DELETE FROM sqlite_sequence WHERE name='clientes'")
    cursor.execute("DELETE FROM sqlite_sequence WHERE name='historico_envios'")
    
    conexao.commit()
    conexao.close()
    print("\n✨ O banco de dados está totalmente zerado e pronto!")

if __name__ == "__main__":
    limpeza_total()