import streamlit as st
import sqlite3
import os

st.set_page_config(page_title="Medical In-Sight | Cadastro", page_icon="🩺")

def salvar_no_banco(nome, email, especialidade, clinica, plano, keywords):
    limites = {"Básico": 1, "Estudante": 2, "Premium": 3}
    tipo_plano = plano.split(" ")[0] 
    limite = limites.get(tipo_plano, 1)
    
    try:
        conexao = sqlite3.connect('medical_insight.db')
        cursor = conexao.cursor()
        # AGORA SALVAMOS TAMBÉM AS KEYWORDS
        cursor.execute('''
            INSERT INTO clientes (nome, email, especialidade, clinica, plano, limite, keywords)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', (nome, email, especialidade, clinica, tipo_plano, limite, keywords))
        conexao.commit()
        conexao.close()
        return True
    except Exception as e:
        st.error(f"Erro ao salvar: {e}")
        return False

st.title("🩺 Medical In-Sight")
st.subheader("Configuração de Curadoria Personalizada")
st.markdown("---")

# --- PARTE DINÂMICA ---
st.write("### 📝 Especialidade e Foco")
lista_especialidades = [
    "Cardiologia", "Psiquiatria", "Pediatria", "Neurologia", 
    "Ginecologia", "Dermatologia", "Ortopedia", "Infectologia", 
    "Endocrinologia", "Oftalmologia", "Oncologia", "Geriatria", 
    "Medicina Interna", "Radiologia", "Anestesiologia", "Outra"
]
escolha = st.selectbox("Selecione sua Especialidade Principal", lista_especialidades)

especialidade_final = escolha
if escolha == "Outra":
    especialidade_final = st.text_input("👉 Digite sua Especialidade:", placeholder="Ex: Nefrologia")

# NOVO CAMPO DE PALAVRAS-CHAVE
keywords = st.text_input("🎯 Temas Específicos (Opcional)", placeholder="Ex: Bipolaridade, TEA, Depressão Resistente")
st.caption("Deixe em branco para receber novidades gerais da sua especialidade.")

# --- FORMULÁRIO ---
with st.form("cadastro_cliente"):
    nome = st.text_input("Nome Completo", placeholder="Ex: Dra. Ana Souza")
    email = st.text_input("E-mail de Recebimento", placeholder="exemplo@gmail.com")
    clinica = st.text_input("Clínica ou Hospital", placeholder="Ex: Hospital das Clínicas")
    
    st.write("---")
    plano = st.radio("Selecione seu Perfil de Atualização:", ["Básico (1 artigo)", "Estudante (2 artigos)", "Premium (3 artigos)"], horizontal=True)
    
    enviar = st.form_submit_button("🚀 Finalizar meu Cadastro")

    if enviar:
        if not nome or not email or (escolha == "Outra" and not especialidade_final):
            st.warning("⚠️ Por favor, preencha os campos obrigatórios (Nome, E-mail e Especialidade).")
        else:
            if salvar_no_banco(nome, email, especialidade_final, clinica, plano, keywords):
                st.success(f"Excelente, {nome}! Cadastro realizado.")
                st.balloons()

if os.path.exists('medical_insight.db'):
    st.sidebar.success("Banco de Dados: Online ✅")