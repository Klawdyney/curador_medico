import streamlit as st
import sqlite3

def atualizar_perfil(email, nova_esp, novas_keys):
    conn = sqlite3.connect('medical_insight.db')
    cursor = conn.cursor()
    # Atualiza especialidade e palavras-chave
    cursor.execute('''
        UPDATE clientes 
        SET especialidade = ?, keywords = ? 
        WHERE email = ?
    ''', (nova_esp, novas_keys, email))
    conn.commit()
    sucesso = cursor.rowcount > 0 
    conn.close()
    return sucesso

def excluir_perfil(email):
    conn = sqlite3.connect('medical_insight.db')
    cursor = conn.cursor()
    cursor.execute('DELETE FROM clientes WHERE email = ?', (email,))
    conn.commit()
    sucesso = cursor.rowcount > 0
    conn.close()
    return sucesso

st.set_page_config(page_title="Medical In-Sight | Gestão", page_icon="🔄")

st.title("🔄 Gestão de Perfil")
st.write("Atualize seus temas de interesse ou gerencie sua conta.")

# --- SEÇÃO DE ATUALIZAÇÃO ---
st.markdown("### 🎯 Ajustar Foco Científico")
email_busca = st.text_input("Confirme seu e-mail cadastrado:")

if email_busca:
    nova_especialidade = st.selectbox("Nova Especialidade Principal:", 
        ["Cardiologia", "Psiquiatria", "Pediatria", "Neurologia", "Nefrologia", "Outra"])
    
    novas_palavras = st.text_input("Novos Temas Específicos (Keywords):", placeholder="Ex: Hipertensão, TEA, Diabetes")
    
    col1, col2 = st.columns([1, 4]) # Alinha o botão à esquerda
    with col1:
        if st.button("Salvar"):
            if atualizar_perfil(email_busca, nova_especialidade, novas_palavras):
                st.success("Perfil atualizado!")
                st.balloons()
            else:
                st.error("E-mail não encontrado.")

    # --- SEÇÃO DE EXCLUSÃO (ÁREA CRÍTICA) ---
    st.markdown("---")
    with st.expander("❌ Encerrar Assinatura e Excluir Dados"):
        st.warning("Atenção: Esta ação é irreversível. Todos os seus dados serão apagados do nosso banco.")
        confirmar_exclusao = st.checkbox("Eu entendo que meus dados serão apagados permanentemente.")
        
        if st.button("Excluir minha conta agora", type="primary", disabled=not confirmar_exclusao):
            if excluir_perfil(email_busca):
                st.success("Cadastro removido com sucesso. Sentiremos sua falta!")
                st.info("Você pode fechar esta aba agora.")
            else:
                st.error("Erro ao excluir. Verifique o e-mail.")