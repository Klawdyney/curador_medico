import streamlit as st
import database_manager as db # Gerente atualizado com UPDATE e busca
from datetime import datetime

# --- 1. CONFIGURAÇÃO DA PÁGINA (DESIGN DE ELITE) ---
st.set_page_config(page_title="Medical In-Sight | Plataforma de Elite", page_icon="🩺", layout="wide")

# --- 2. LISTAS E DICIONÁRIOS (SINCRO TOTAL SITE/ROBÔ) ---
ESPECIALIDADES_OFICIAIS = [
    "Cardiologia", "Psiquiatria", "Pediatria", "Neurologia", "Dermatologia", 
    "Oftalmologia", "Cirurgia Geral", "Ginecologia e Obstetrícia", "Radiologia",
    "Anestesiologia", "Oncologia", "Endocrinologia", "Ortopedia", "Medicina Interna", "Outra"
]

LIMITES_PLANOS = {
    "Básico": 1,
    "Estudante": 2,
    "Premium": 5,
    "Sob Demanda": 10
}

DIAS_NOMES = {
    "Segunda-feira": "seg", "Terça-feira": "ter", "Quarta-feira": "qua",
    "Quinta-feira": "qui", "Sexta-feira": "sex", "Sábado": "sab", "Domingo": "dom"
}

# --- 3. ESTILO VISUAL CUSTOMIZADO (CSS PROFISSIONAL) ---
st.markdown("""
    <style>
    .main { background-color: #f8f9fa; }
    [data-testid="stSidebar"] button[data-baseweb="tab"] div p {
        font-size: 20px !important;
        font-weight: bold !important;
    }
    .stButton>button { 
        width: 100%; border-radius: 8px; height: 3.5em; 
        background-color: #005ea2; color: white; font-weight: 600; border: none; 
    }
    .banner-container {
        width: 100%; height: 280px; overflow: hidden; border-radius: 12px;
        margin-bottom: 20px; box-shadow: 0 4px 10px rgba(0,0,0,0.1);
    }
    .banner-container img {
        width: 100%; height: 100%; object-fit: cover;
    }
    .plan-card {
        background-color: white; padding: 25px; border-radius: 15px;
        border: 1px solid #e0e0e0; text-align: center; height: 100%;
    }
    </style>
    """, unsafe_allow_html=True)

# --- 4. CENTRAL DE TRANSPARÊNCIA (FAQ COMPLETO) ---
def exibir_faq_completo():
    st.divider()
    st.subheader("🛡️ Central de Transparência e Rigor Científico")
    
    with st.expander("🔍 1. VERACIDADE E QUALIDADE DAS INFORMAÇÕES"):
        st.write("Nosso sistema conecta-se diretamente à API do PubMed/MEDLINE. Não utilizamos blogs; apenas artigos com DOI rastreável.")

    with st.expander("🎯 2. OBJETIVO E FINALIDADE DO SERVIÇO"):
        st.write("Nossa finalidade é a curadoria e otimização de tempo, filtrando o ruído para que você foque no que realmente impacta sua conduta clínica.")

    with st.expander("🔐 3. SEGURANÇA E PRIVACIDADE (LGPD)"):
        st.write("Operamos sob os pilares da LGPD. Seus dados e keywords são criptografados.")

    with st.expander("📂 4. GESTÃO DA ASSINATURA E FIDELIDADE"):
        st.write("No Portal do Assinante, você tem autonomia total para atualizar seus focos de interesse. Não há fidelidade.")

    with st.expander("📊 5. QUALIDADE DO RESUMO E DADOS"):
        st.write("Nossa IA utiliza a técnica 'Chain-of-Thought', extraindo Metodologia, Resultados e P-valores.")

    with st.expander("🌐 6. FONTES E DIRETRIZES"):
        st.write("Focamos no PubMed (padrão ouro), mas também identificamos diretrizes de sociedades como AHA, ESC e AAD.")

    with st.expander("🛠️ 7. SUPORTE TÉCNICO E PERSONALIZAÇÃO"):
        st.write("Suporte via e-mail e WhatsApp em horário comercial para ajustes finos em suas keywords.")

    with st.expander("🤖 8. TRATAMENTO DE ALUCINAÇÕES DE IA"):
        st.write("Sempre fornecemos o link direto (DOI) da fonte original em cada resumo para conferência.")

# --- 5. NAVEGAÇÃO LATERAL ---
aba_home, aba_login, aba_assinar = st.sidebar.tabs(["Home", "Login", "Assinar"])

# --- 6. ABA DE ASSINATURA (COM KEYWORDS OPCIONAIS) ---
with aba_assinar:
    st.subheader("Nova Assinatura Profissional")
    with st.form("registro_vendas", clear_on_submit=False):
        nome_n = st.text_input("Nome Completo:")
        email_n = st.text_input("E-mail Profissional (Login):")
        whatsapp_n = st.text_input("WhatsApp:", placeholder="5531999999999", help="55 + DDD + Número")
        
        esp_n = st.selectbox("Especialidade Principal:", ESPECIALIDADES_OFICIAIS)
        
        # Keywords agora é opcional no momento do cadastro
        keywords_n = st.text_input("Foco / Keywords (Opcional):", placeholder="Ex: Melanoma, Glaucoma...")
        
        plano_n = st.selectbox("Selecione seu Plano:", list(LIMITES_PLANOS.keys()))
        
        st.write("--- Configurações de Recebimento ---")
        dia_n_bonito = st.selectbox("Melhor dia para envio:", list(DIAS_NOMES.keys()))
        dia_n = DIAS_NOMES[dia_n_bonito]
        hora_n = st.selectbox("Melhor horário:", ["08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00", "16:00", "17:00", "18:00", "19:00"])
        
        if st.form_submit_button("Confirmar Assinatura e Ativar IA"):
            if nome_n and email_n and whatsapp_n:
                kw_final = keywords_n if keywords_n else esp_n # Padrão é a especialidade
                limite_artigos = LIMITES_PLANOS.get(plano_n, 1)
                
                sucesso = db.cadastrar_medico({
                    "nome": nome_n, "email": email_n, "whatsapp": whatsapp_n,
                    "especialidade": esp_n, "clinica": "Particular", "keywords": kw_final,
                    "plano": plano_n, "limite": limite_artigos, 
                    "dia_envio": dia_n, "horario_envio": hora_n
                })
                if sucesso:
                    st.success(f"✅ Bem-vindo, Dr. {nome_n}! Seu plano {plano_n} está ativo.")
                    st.balloons()
            else:
                st.error("Por favor, preencha todos os campos obrigatórios.")

# --- 7. ABA DE LOGIN E PAINEL EDITÁVEL (A VIRADA DE CHAVE) ---
with aba_login:
    st.write("Acesso restrito para assinantes.")
    email_login = st.text_input("E-mail cadastrado:", key="login_field")

if email_login:
    perfil = db.buscar_medico_por_email(email_login) # Busca dados reais no Supabase
    
    if perfil:
        st.success(f"Painel do Assinante: Conectado como {email_login}")
        tab_perfil, tab_biblioteca, tab_ajuda = st.tabs(["⚙️ Meu Perfil", "📚 Biblioteca de PDFs", "🛡️ Ajuda"])
        
        with tab_perfil:
            st.subheader("Personalize sua Inteligência Curadora")
            with st.form("update_perfil"):
                c1, c2 = st.columns(2)
                with c1:
                    novo_plano = st.selectbox("Mudar Plano:", list(LIMITES_PLANOS.keys()), 
                                              index=list(LIMITES_PLANOS.keys()).index(perfil['plano']))
                    nova_esp = st.selectbox("Atualizar Especialidade:", ESPECIALIDADES_OFICIAIS,
                                            index=ESPECIALIDADES_OFICIAIS.index(perfil['especialidade']))
                with c2:
                    dia_sigla = perfil['dia_envio']
                    dia_nome_atual = [k for k, v in DIAS_NOMES.items() if v == dia_sigla][0]
                    novo_dia_n = st.selectbox("Dia de Recebimento:", list(DIAS_NOMES.keys()), 
                                              index=list(DIAS_NOMES.keys()).index(dia_nome_atual))
                    novo_h = st.selectbox("Horário de Recebimento:", ["08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00", "16:00", "17:00", "18:00", "19:00"],
                                          index=["08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00", "16:00", "17:00", "18:00", "19:00"].index(perfil['horario_envio']))
                
                novas_kw = st.text_area("Refinar Foco de Busca (Keywords):", value=perfil['keywords'], help="Separe termos por vírgula para refinar a IA.")
                
                if st.form_submit_button("Salvar Alterações no Perfil"):
                    sucesso_up = db.atualizar_perfil_medico(email_login, {
                        "plano": novo_plano,
                        "limite": LIMITES_PLANOS[novo_plano],
                        "especialidade": nova_esp,
                        "keywords": novas_kw,
                        "dia_envio": DIAS_NOMES[novo_dia_n],
                        "horario_envio": novo_h,
                        "clinica": perfil['clinica']
                    })
                    if sucesso_up:
                        st.success("✅ Perfil atualizado! As mudanças já valem para o próximo ciclo do robô.")
                        st.rerun()
    else:
        st.warning("E-mail não identificado em nossa base ativa.")

    with tab_biblioteca:
        st.info("Aqui aparecerão seus boletins semanais em PDF para download.")
    with tab_ajuda: exibir_faq_completo()

else:
    # --- 8. HOME PAGE (ESTRUTURA ORIGINAL INTEGRAL) ---
    st.markdown(f"""
        <div class="banner-container">
            <img src="https://images.unsplash.com/photo-1576091160550-2173dba999ef?q=80&w=2070">
        </div>
    """, unsafe_allow_html=True)
    
    c_msg, c_stats = st.columns([2, 1])
    with c_msg:
        st.header("Inteligência Científica, otimize seu tempo clínico.")
        st.write("Curadoria automatizada do PubMed entregue no seu E-mail e WhatsApp.")
        st.link_button("📄 Ver Exemplo de Boletim", "https://seu-exemplo.pdf")
    
    with c_stats:
        st.metric("Base PubMed", "+3.000 Artigos/Dia")
        st.metric("Rigor IA", "Protocolo Grounded AI")

    st.divider()
    st.subheader("Nossos Planos")
    p1, p2, p3, p4 = st.columns(4)
    p1.markdown("<div class='plan-card'><b>BÁSICO</b><br>R$ 39,90/mês<br><small>1 Artigos Semanais</small></div>", unsafe_allow_html=True)
    p2.markdown("<div class='plan-card'><b>ESTUDANTE</b><br>R$ 49,90/mês<br><small>2 Artigos Semanal</small></div>", unsafe_allow_html=True)
    p3.markdown("<div class='plan-card'><b>PREMIUM</b><br>R$ 79,90/mês<br><small>5 Artigos Semanais</small></div>", unsafe_allow_html=True)
    p4.markdown("<div class='plan-card'><b>SOB DEMANDA</b><br>Consulte<br><small>Clínicas/Grupos</small></div>", unsafe_allow_html=True)

    exibir_faq_completo()

st.divider()
st.caption(f"Medical In-Sight © 2026 | Tecnologia em ADS | {datetime.now().year}")