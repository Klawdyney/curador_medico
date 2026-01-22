import time
from google import genai
from Bio import Entrez
from fpdf import FPDF

# --- CONFIGURAÇÕES ---
MINHA_CHAVE = "AIzaSyBAIyjlLexxUtFazc7NxSXqbl5fmFCYVAE"
Entrez.email = "claudinei.jb@gmail.com"
client = genai.Client(api_key=MINHA_CHAVE)

class PDF_Premium(FPDF):
    def header(self):
        self.set_fill_color(0, 51, 102) 
        self.rect(0, 0, 210, 35, 'F')
        self.set_font("helvetica", 'B', 22)
        self.set_text_color(255, 255, 255)
        self.cell(0, 15, "IN-SIGHT: Boletim Cientifico", new_x="LMARGIN", new_y="NEXT", align='C')
        self.set_font("helvetica", 'I', 10)
        self.cell(0, 5, "Curadoria Especializada de Alta Performance", new_x="LMARGIN", new_y="NEXT", align='C')
        self.ln(15)

    def footer(self):
        self.set_y(-15)
        self.set_font("helvetica", 'I', 8)
        self.set_text_color(128, 128, 128)
        self.cell(0, 10, f"Pagina {self.page_no()} | Fonte: PubMed | IA: Gemini 2.0 Flash", align='C')

def buscar_multiplos_artigos(tema, quantidade=3):
    print(f"\n[1/3] Consultando base PubMed para: {tema}...")
    try:
        handle = Entrez.esearch(db="pubmed", term=tema, retmax=quantidade, sort="date")
        record = Entrez.read(handle)
        handle.close()
        artigos = []
        if record["IdList"]:
            for id_artigo in record["IdList"]:
                fetch = Entrez.efetch(db="pubmed", id=id_artigo, rettype="abstract", retmode="text")
                artigos.append(fetch.read())
                fetch.close()
            return artigos
    except Exception as e:
        print(f"Erro na busca: {e}")
    return []

def gerar_relatorio_visual(conteudo, especialidade):
    pdf = PDF_Premium()
    pdf.set_auto_page_break(auto=True, margin=20)
    pdf.add_page()
    pdf.ln(10)
    
    data_atual = time.strftime("%d/%m/%Y")
    pdf.set_font("helvetica", 'B', 11)
    pdf.set_text_color(0, 51, 102)
    pdf.cell(0, 10, f"EDICAO: {data_atual}   |   ESPECIALIDADE: {especialidade.upper()}", 
             new_x="LMARGIN", new_y="NEXT", align='L')
    pdf.ln(5)
    
    pdf.set_font("helvetica", size=11)
    pdf.set_text_color(0, 0, 0)
    texto_seguro = conteudo.encode('latin-1', 'replace').decode('latin-1')
    pdf.multi_cell(0, 8, txt=texto_seguro)
    
    nome_arquivo = f"Boletim_{especialidade}_{int(time.time())}.pdf"
    pdf.output(nome_arquivo)
    return nome_arquivo

def exibir_menu():
    print("\n" + "="*30)
    print("  MEDICAL IN-SIGHT - MENU  ")
    print("="*30)
    print("1. Psiquiatria")
    print("2. Neurologia")
    print("3. Pediatria")
    print("4. Cardiologia")
    print("5. Digitar outro tema...")
    print("0. Sair")
    return input("\nEscolha a area medica: ")

def main():
    while True:
        opcao = exibir_menu()
        
        if opcao == '0':
            print("Encerrando sistema. Ate amanha!")
            break
        
        temas = {
            '1': "Psychiatry latest clinical trials 2026",
            '2': "Neurology breakthrough research 2026",
            '3': "Pediatrics clinical updates 2026",
            '4': "Cardiology new guidelines 2026"
        }
        
        if opcao in temas:
            tema_busca = temas[opcao]
            nome_especialidade = temas[opcao].split()[0] # Pega a primeira palavra
        elif opcao == '5':
            nome_especialidade = input("Digite a especialidade (ex: Dermatologia): ")
            tema_busca = f"{nome_especialidade} latest research 2026"
        else:
            print("Opcao invalida.")
            continue

        lista = buscar_multiplos_artigos(tema_busca, quantidade=3)
        
        if lista:
            print(f"[2/3] Gemini 2.0 analisando evidencias...")
            texto_unificado = "\n\n--- ESTUDO ---\n".join(lista)
            prompt = f"""Analise estes 3 artigos de {nome_especialidade}. 
            Crie um boletim executivo em Portugues com:
            - Titulo Traduzido
            - Resumo Clinico (O que muda na pratica)
            - Conclusao
            Estudos: {texto_unificado[:8000]}"""
            
            response = client.models.generate_content(model="gemini-2.0-flash", contents=prompt)
            arquivo = gerar_relatorio_visual(response.text, nome_especialidade)
            print(f"\n[OK] Boletim gerado: {arquivo}")
        else:
            print("Nenhum artigo encontrado.")

if __name__ == "__main__":
    main()