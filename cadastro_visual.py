import streamlit as st
import sqlite3

# 1. INICIALIZAÇÃO DA MEMÓRIA (Session State)
# Isso garante que os dados fiquem salvos mesmo se a página recarregar por causa de um erro
if "nome_temp" not in st.session_state: st.session_state.nome_temp = ""
if "email_temp" not in st.session_state: st.session_state.email_temp = ""
if "whatsapp_temp" not in st.session_state: st.session_state.whatsapp_temp = "55"

# Configuração da Página
st.set_page_config(page_title="Cadastro Oficial | Medical In-Sight", page_icon="🩺", layout="centered")

st.title("🩺 Cadastro de Novo Assinante")
st.markdown("""
    Use esta interface para cadastrar médicos. 
    **Nota:** Os campos corretos serão preservados caso haja algum erro de preenchimento.
""")

# 2. FORMULÁRIO (clear_on_submit deve ser FALSE para controlarmos manualmente)
with st.form("form_oficial", clear_on_submit=False):
    st.subheader("Informações Pessoais")
    
    # Vinculamos o conteúdo dos inputs às variáveis da memória (session_state)
    nome = st.text_input("Nome Completo:", value=st.session_state.nome_temp)
    email = st.text_input("E-mail de Acesso:", value=st.session_state.email_temp)
    
    whatsapp = st.text_input(
        "WhatsApp (DDD + Número):", 
        value=st.session_state.whatsapp_temp, 
        max_chars=13,
        help="Digite apenas números: 55 + DDD + Número. Exemplo: 5511999998888"
    )

    st.markdown("---")
    st.subheader("Configurações de Curadoria")
    
    especialidade = st.selectbox("Especialidade:", 
        ["Cardiologia", "Psiquiatria", "Pediatria", "Neurologia", "Nefrologia", 
        "Ginecologia e Obstetrícia", "Dermatologia", "Ortopedia", "Endocrinologia", 
        "Oncologia", "Geriatria", "Medicina Interna", "Medicina de Família", 
        "Infectologia", "Gastrenterologia", "Urologia", "Oftalmologia", "Outra"])
    
    clinica = st.text_input("Clínica/Hospital:", value="Consultório Particular")
    
    keywords = st.text_area("Keywords de Interesse (separe por vírgula):", 
                            placeholder="Ex: Heart Failure, Guidelines, Immunotherapy")
    
    col1, col2 = st.columns(2)
    with col1:
        plano = st.selectbox("Plano Contratado:", ["Básico", "Premium", "Estudante"])
    with col2:
        dia_envio = st.selectbox("Dia do Boletim:", ["Segunda", "Terça", "Quarta", "Quinta", "Sexta", "Sábado", "Domingo"])

    submit = st.form_submit_button("Finalizar Cadastro")

    if submit:
        # --- 3. VALIDAÇÃO E LIMPEZA SELETIVA ---
        whatsapp_limpo = "".join(filter(str.isdigit, whatsapp))
        tem_erro = False

        # Validação do Nome
        if len(nome.strip()) < 3:
            st.error("⚠️ Nome inválido ou muito curto. O campo foi limpo.")
            st.session_state.nome_temp = "" # Apaga apenas o nome na memória
            tem_erro = True
        else:
            st.session_state.nome_temp = nome # Mantém o que está certo

        # Validação do E-mail
        if "@" not in email or "." not in email:
            st.error(f"⚠️ O e-mail '{email}' é inválido. O campo foi limpo.")
            st.session_state.email_temp = "" # Apaga apenas o e-mail
            tem_erro = True
        else:
            st.session_state.email_temp = email

        # Validação do WhatsApp
        if not whatsapp_limpo.startswith("55") or len(whatsapp_limpo) < 12:
            st.error("⚠️ WhatsApp incorreto. Use: 55 + DDD + Número. O campo foi resetado.")
            st.session_state.whatsapp_temp = "55" # Reseta apenas o WhatsApp
            tem_erro = True
        else:
            st.session_state.whatsapp_temp = whatsapp_limpo

        # Se encontrou erro em qualquer campo, para e recarrega a tela com os campos limpos
        if tem_erro:
            st.rerun() 
        
        else:
            # --- 4. TUDO OK, SALVAR NO BANCO ---
            try:
                if plano == "Premium": limite_artigos = 3
                elif plano == "Estudante": limite_artigos = 2
                else: limite_artigos = 1

                conn = sqlite3.connect('medical_insight.db')
                cursor = conn.cursor()
                cursor.execute("""
                    INSERT INTO clientes (
                        nome, email, whatsapp, especialidade, clinica, 
                        keywords, plano, limite, dia_envio, horario_envio
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (nome, email, whatsapp_limpo, especialidade, clinica, 
                      keywords, plano, limite_artigos, dia_envio, "08:00"))
                
                conn.commit()
                conn.close()
                
                # LIMPEZA TOTAL APÓS SUCESSO
                st.session_state.nome_temp = ""
                st.session_state.email_temp = ""
                st.session_state.whatsapp_temp = "55"
                
                st.success(f"🚀 Dr(a). {nome} cadastrado com sucesso!")
                st.balloons()
                st.rerun() # Atualiza a tela para mostrar os campos vazios após o sucesso

            except sqlite3.IntegrityError:
                st.error("❌ Este e-mail já está cadastrado!")
            except Exception as e:
                st.error(f"⚠️ Erro ao salvar: {e}")

# Rodapé
st.markdown("---")
st.caption("Medical In-Sight Admin | Inteligência de Validação de Campos v2.0")