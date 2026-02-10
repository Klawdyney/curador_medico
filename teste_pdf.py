import logging
from app import processar_medico_completo
from dotenv import load_dotenv

# Carrega as configurações do .env
load_dotenv()

# Dados para o teste de elite
meu_usuario_teste = {
    "nome": "Claudinei Barros",
    "email": "claudinei.jb@gmail.com", # O PDF será enviado para aqui se o SMTP estiver ok
    "especialidade": "Psiquiatria",
    "keywords": "Depressão, Tratamento, Psicopatologia",
    "whatsapp": "5531994007459",
    "clinica": "Medical In-Sight Premium",
    "limite": 2
}

if __name__ == "__main__":
    print("🚀 Acionando o motor Premium (sem Streamlit)...")
    resultado = processar_medico_completo(meu_usuario_teste)
    print(f"🏁 Resultado do Processamento: {resultado}")
    print("\n✅ Verifique a pasta do projeto: Boletim_Claudinei_Barros.pdf")