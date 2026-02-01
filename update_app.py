import streamlit as st
import sqlite3

def criar_conexao():
    return sqlite3.connect('medical_insight.db')

# --- FUNÇÕES DE BANCO DE DADOS ATUALIZADAS ---
def buscar_dados_atuais(email):
    conn = criar_conexao()
    cursor = conn.cursor()
    # Adicionamos os novos campos na busca
    cursor.execute('''
        SELECT nome, especialidade, keywords, whatsapp, plano, limite, dia_envio, horario_envio, valor_assinatura 
        FROM clientes WHERE email = ?
    ''', (email,))
    resultado = cursor.fetchone()
    conn.close()
    return resultado

def atualizar_perfil(email, nova_esp, novas_keys, novo_zap, novo_plano, novo_limite, novo_dia, novo_horario, novo_valor):
    conn = criar_conexao()
    cursor = conn.cursor()
    cursor.execute('''
        UPDATE clientes 
        SET especialidade = ?, keywords = ?, whatsapp = ?, plano = ?, limite = ?, dia_envio = ?, horario_envio = ?, valor_assinatura = ?
        WHERE email = ?
    ''', (nova_esp, novas_keys, novo_zap, novo_plano, novo_limite, novo_dia, novo_horario, novo_valor, email))
    conn.commit()
    sucesso = cursor.rowcount > 0 
    conn.close()
    return sucesso

def cadastrar_novo_medico(nome, email, esp, keys, zap, plano, limite, dia, horario, valor):
    conn = criar_conexao()
    cursor = conn.cursor()
    try:
        cursor.execute('''
            INSERT INTO clientes (nome, email, especialidade, keywords, whatsapp, plano, limite, dia_envio, horario_envio, valor_assinatura)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (nome, email, esp, keys, zap, plano, limite, dia, horario, valor))
        conn.commit()
        return True
    except Exception as e:
        st.error(f"Erro ao cadastrar: {e}")
        return False
    finally:
        conn.close()

def excluir_perfil(email):
    conn = criar_conexao()
    cursor = conn.cursor()
    cursor.execute('DELETE FROM clientes WHERE email = ?', (email,))
    conn.commit()
    sucesso = cursor.rowcount > 0
    conn.close()
    return sucesso

# --- CONFIGURAÇÃO DA INTERFACE ---
st.set_page_config(page_title="Medical In-Sight | Gestão", page_icon="🩺")
st.title("🩺 Painel de Controle Medical In-Sight")

aba1, aba2 = st.tabs(["📝 Novo Cadastro", "🔄 Atualizar/Excluir"])

lista_especialidades = [
    "Cardiologia", "Psiquiatria", "Pediatria", "Neurologia", "Nefrologia", 
    "Ginecologia e Obstetrícia", "Dermatologia", "Ortopedia", "Endocrinologia", 
    "Oncologia", "Geriatria", "Medicina Interna", "Medicina de Família", 
    "Infectologia", "Gastrenterologia", "Urologia", "Oftalmologia", "Outra"
]

lista_dias = ["Segunda", "Terça", "Quarta", "Quinta", "Sexta", "Sábado", "Domingo", "Diário"]

# --- ABA 1: NOVO CADASTRO ---
with aba1:
    st.subheader("Cadastrar Novo Médico")
    with st.form("form_cadastro"):
        col1, col2 = st.columns(2)
        with col1:
            nome_novo = st.text_input("Nome do Médico:")
            email_novo = st.text_input("E-mail:")
            esp_novo = st.selectbox("Especialidade:", lista_especialidades)
        with col2:
            keys_novo = st.text_input("Temas (Keywords):", placeholder="Ex: Hipertensão, TEA")
            zap_novo = st.text_input("WhatsApp:", value="55", help="Formato: 55 + DDD + Número")
            horario_novo = st.time_input("Horário de Entrega:", value=None)
        
        st.divider()
        st.write("**Selecione o Plano de Assinatura:**")
        # Novos planos conforme sua decisão
        plano_op = st.radio(
            "Planos disponíveis:", 
            ["Básico (R$ 89,90 | 1x semana | 2 artigos.)", 
             "Intermediário (R$ 159,90 | 2x semana | 4 artigos.)", 
             "Elite (R$ 249,90 | Diário | 6 artigos.)"], 
            horizontal=True
        )
        
        mapa_planos = {
            "Básico (R$ 89,90 | 1x semana | 2 artigos.)": ("Básico", 2, 89.90),
            "Intermediário (R$ 159,90 | 2x semana | 4 artigos.)": ("Intermediário", 4, 159.90),
            "Elite (R$ 249,90 | Diário | 6 artigos.)": ("Elite", 6, 249.90)
        }

        dia_novo = st.selectbox("Dia de preferência para envio:", lista_dias)
        
        if st.form_submit_button("Finalizar Cadastro"):
            zap_limpo = ''.join(filter(str.isdigit, zap_novo))
            email_valido = "@" in email_novo and "." in email_novo
            
            if not nome_novo or not email_novo or not horario_novo:
                st.error("⚠️ Nome, E-mail e Horário são obrigatórios.")
            elif not email_valido:
                st.error("⚠️ E-mail inválido.")
            else:
                p_nome, p_limite, p_valor = mapa_planos[plano_op]
                h_formatado = horario_novo.strftime("%H:%M")
                
                if cadastrar_novo_medico(nome_novo, email_novo, esp_novo, keys_novo, zap_limpo, p_nome, p_limite, dia_novo, h_formatado, p_valor):
                    st.success(f"Dr(a). {nome_novo} cadastrado com sucesso no Plano {p_nome}!")
                    st.balloons()

# --- ABA 2: ATUALIZAÇÃO (RESUMIDA PARA O NOVO MODELO) ---
with aba2:
    st.subheader("Gerenciar Médico Existente")
    email_busca = st.text_input("Digite o e-mail cadastrado para buscar:")

    if email_busca:
        dados = buscar_dados_atuais(email_busca)
        if dados:
            nome_at, esp_at, keys_at, zap_at, plano_at, limite_at, dia_at, hora_at, valor_at = dados
            st.info(f"Editando perfil de: **{nome_at}** | Plano Atual: **{plano_at}**")
            
            col1, col2 = st.columns(2)
            with col1:
                nova_esp = st.selectbox("Especialidade:", lista_especialidades, index=lista_especialidades.index(esp_at) if esp_at in lista_especialidades else 0)
                novas_keys = st.text_input("Keywords:", value=keys_at)
                novo_zap = st.text_input("WhatsApp:", value=zap_at)
            with col2:
                novo_dia = st.selectbox("Dia de Envio:", lista_dias, index=lista_dias.index(dia_at) if dia_at in lista_dias else 0)
                novo_horario = st.text_input("Horário (HH:MM):", value=hora_at)
            
            st.divider()
            novo_plano_escolha = st.radio("Alterar Plano:", ["Básico", "Intermediário", "Elite"], 
                                          index=["Básico", "Intermediário", "Elite"].index(plano_at) if plano_at in ["Básico", "Intermediário", "Elite"] else 0)
            
            mapa_up = {"Básico": (2, 89.90), "Intermediário": (4, 159.90), "Elite": (6, 249.90)}

            if st.button("Salvar Alterações"):
                limite_up, valor_up = mapa_up[novo_plano_escolha]
                if atualizar_perfil(email_busca, nova_esp, novas_keys, novo_zap, novo_plano_escolha, limite_up, novo_dia, novo_horario, valor_up):
                    st.success("Dados atualizados com sucesso!")
                    st.rerun()
        else:
            if email_busca: st.warning("Médico não encontrado.")

    st.markdown("---")
    with st.expander("❌ Área Crítica: Excluir Cadastro"):
        if st.button("Excluir conta agora", type="primary"):
            if excluir_perfil(email_busca):
                st.success("Cadastro removido.")
                st.rerun()