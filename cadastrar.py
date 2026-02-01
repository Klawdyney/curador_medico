import sqlite3

def novo_cadastro():
    print("\n--- CADASTRO DE NOVO CLIENTE MEDICAL IN-SIGHT ---")
    nome = input("Nome Completo: ")
    email = input("E-mail: ")
    
    # LISTA DE ESPECIALIDADES PADRONIZADA
    especialidades = [
        "Cardiologia", "Psiquiatria", "Pediatria", "Neurologia", 
        "Nefrologia", "Ginecologia e Obstetrícia", "Dermatologia", 
        "Ortopedia", "Endocrinologia", "Oncologia", "Geriatria"
    ]
    
    print("\nEspecialidades disponíveis:")
    for i, esp in enumerate(especialidades, 1):
        print(f"{i}. {esp}")
    
    op_esp = input("Escolha o número da especialidade (ou digite o nome de outra): ")
    
    # Se o usuário digitou um número, pega da lista. Se digitou texto, usa o texto.
    if op_esp.isdigit() and 1 <= int(op_esp) <= len(especialidades):
        especialidade = especialidades[int(op_esp)-1]
    else:
        especialidade = op_esp

    clinica = input("Nome da Clínica/Hospital: ")
    print("\nPlanos: 1. Basico (1 artigo) | 2. Premium (3 artigos) | 3. Estudante (2 artigos)")
    opcao_plano = input("Escolha o Plano (1, 2 ou 3): ")
    
    planos = {"1": ("Basico", 1), "2": ("Premium", 3), "3": ("Estudante", 2)}
    plano_nome, limite = planos.get(opcao_plano, ("Basico", 1))

    try:
        conexao = sqlite3.connect('medical_insight.db')
        cursor = conexao.cursor()
        cursor.execute('''
            INSERT INTO clientes (nome, email, especialidade, clinica, plano, limite)
            VALUES (?, ?, ?, ?, ?, ?)
        ''', (nome, email, especialidade, clinica, plano_nome, limite))
        conexao.commit()
        conexao.close()
        print(f"\n[SUCESSO] {nome} cadastrado com sucesso!")
    except Exception as e:
        print(f"\n[ERRO] Não foi possível cadastrar: {e}")

if __name__ == "__main__":
    while True:
        novo_cadastro()
        continuar = input("\nDeseja cadastrar outro? (s/n): ")
        if continuar.lower() != 's':
            break