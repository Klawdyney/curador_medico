import sqlite3

def deletar_usuario(nome_aproximado):
    conexao = sqlite3.connect('medical_insight.db')
    cursor = conexao.cursor()
    
    # Deleta qualquer usuário que tenha esse nome
    cursor.execute("DELETE FROM clientes WHERE nome LIKE ?", (f'%{nome_aproximado}%',))
    
    conexao.commit()
    print(f"Registros contendo '{nome_aproximado}' foram removidos!")
    conexao.close()

if __name__ == "__main__":
    # Digite o nome ou parte do nome que quer apagar
    nome = input("Digite o nome do cliente para apagar: ")
    deletar_usuario(nome)