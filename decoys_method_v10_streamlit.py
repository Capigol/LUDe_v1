# -*- coding: utf-8 -*-
"""
Created on Mon May 10 14:27:11 2021

@author: Lucas
"""

##### Decoys generator #####
### This script identify and select molecules into ChEMBL27 database 
### these are physicochemical similar and also topological/structural different to the uploaded compounds 


# Needed packages
from pathlib import Path
import streamlit as st
import pandas as pd
import base64
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs # para calcular Tanimoto similarity
from rdkit.Chem.Scaffolds import MurckoScaffold # para poder calcular los frameworks
from rdkit.Chem import Descriptors # para poder calcular descriptores
from rdkit.Chem import rdFMCS # para calcular la MCS
from molvs import standardize_smiles
import random
from time import process_time
import os

#---------------------------------#
# Page layout
## Page expands to full width
st.set_page_config(page_title='LIDEB Tools - Decoys',
    layout='wide')

######
# Function to put a picture as header   
def img_to_bytes(img_path):
    img_bytes = Path(img_path).read_bytes()
    encoded = base64.b64encode(img_bytes).decode()
    return encoded

from PIL import Image
image = Image.open('cropped-header-lude.png')
st.image(image)

#header_html = "<img src='data:image/png;base64,{}' class='img-fluid'>".format(img_to_bytes("cropped-header-lude.png"))
#st.markdown(header_html, unsafe_allow_html=True,)
######

st.markdown("![Twitter Follow](https://img.shields.io/twitter/follow/LIDeB_UNLP?style=social)")
st.subheader(":pushpin:" "About Us")
st.markdown("We are a team interested to develop new tools cheminformatics for using in areas of computer-assisted drug design and machine learning in drug discovery. We belong to Laboratory of Bioactive Research and Development (LIDeB) of the National University of La Plata, Argentina. Our research group is focused on rational designing of new drug candidates and drug repurposing for epilepsy and neglected tropical diseases such as Chagas disease, leishmaniasis, malaria.Another important goal of our group is the development and caracterization of nanocarriers. The work developed by our group has resulted in publications in international indexed journals, abstracts, congress and awards in national and international of Medicinal and Computational Chemistry.")
st.markdown(":computer:""**Web Site** " "<https://lideb.biol.unlp.edu.ar>")

# Start the stopwatch / counter 
t1_start = process_time() 


#---------------------------------#
st.write("""
# LIDeB Tools - DECOYS WEB APP

**WebApp to generate decoys from a set of active compounds**

**Decoys are molecules that have not been tested against a molecular target but due to its structural features it should not bind to the target.**
 They are very usefull for validation of machine learning models.

In this Web App, decoys are obtained by three sequential steps:
    
- Searching molecules with similar physicochemical properties of the loaded molecules in a curated ChEMBL database.
- Filtering the selected molecules by dissimilarity against each individual loaded molecule.
- Randomly selecting a desired number of decoys for each individual loaded molecule.
- Filtering the selected molecules by the dissimilarity against all the loaded molecules.
Finally, you can download a file with your **decoys.**

The next workflow summarizes the steps performed by this method:
""")

text = '''
---

'''


#st.markdown("""
 #        **To cite the application, please reference XXXXXXXXX**
 #        """)


image = Image.open('workflow_lude.png')
st.image(image, caption='LUDe Workflow')

st.markdown(text)


# OPTIONS

#---------------------------------#
# Sidebar - Collects user input features into dataframe
st.sidebar.header('Upload your smiles')

uploaded_file_1 = st.sidebar.file_uploader("Upload your smiles in a TXT file", type=["txt"])

st.sidebar.markdown("""
[Example TXT input file](moleculew.txt)
""")

lude_setting = st.sidebar.checkbox('Check to change the default configuration', value=False)
if lude_setting == True:
    # Sidebar - Physicochemical features
    st.sidebar.header('Physicochemical features limits')
    lim_MW = st.sidebar.slider('Molecular weight', 0, 10, 5, 1)
    lim_logP = st.sidebar.slider('LogP', 0.0, 0.5, 0.5, 0.1)
    lim_rot_bonds = st.sidebar.slider('Rotable bonds', 0, 3, 1, 1)
    lim_N_Haccep = st.sidebar.slider('Num of H Acceptors', 0, 2, 1, 1)
    lim_N_Hdon = st.sidebar.slider('Num of H Donors', 0, 2, 1, 1)
    
    # Sidebar - Dissimilarty conditions
    st.sidebar.subheader('Dissimilarty conditions')
    
    fingerprint_radio = st.sidebar.slider('Set the fingerprint radio', 1, 3, 2, 1)
    fingerprint_lenght = st.sidebar.slider('Set the fingerprint lenght', 512, 2048, 1024, 512)
    similarity_metric = st.sidebar.selectbox("Select the similarity metric", ("TanimotoSimilarity", "DiceSimilarity", "CosineSimilarity", "SokalSimilarity", "RusselSimilarity", "KulczynskiSimilarity", "McConnaugheySimilarity"),0)
    
    max_similarity_limit = st.sidebar.slider('Maximum allowed similarity value', 0.0, 1.0, 0.2, 0.05)
    
    # Limit of the fraction of the Maximum Common Substructure
    lim_fraction_MCS = st.sidebar.slider('Limit of the fraction of the Maximum Common Substructure', 0.0, 1.0, 0.5, 0.05)
    
    # Decoys with different framework. Options "Yes" or "No"
    framework_filter = st.sidebar.checkbox('Framework filter',value=True)
    
    # Maximum similarity allowed between decoys and any of the actives
    max_similarity_limit_all = st.sidebar.slider('Maximum similarity allowed between decoys and any of the actives', 0.0, 0.7, 0.3, 0.05)
    
    # Sidebar - Max number of decoys
    st.sidebar.subheader('Maximum number of decoys')
    
    # Max number of decoys by active
    max_decoys = st.sidebar.slider('Maximum number of decoys by active compound', 50, 1000, 100, 50)

else:
    lim_MW = 5
    lim_logP = 0.5
    lim_rot_bonds = 1
    lim_N_Haccep = 1
    lim_N_Hdon = 1
    fingerprint_radio = 2
    fingerprint_lenght = 1024
    similarity_metric ="TanimotoSimilarity"
    max_similarity_limit = 0.2
    lim_fraction_MCS = 0.5
    framework_filter = True
    max_similarity_limit_all = 0.3
    max_decoys = 100
    
#---------------------------------#
# Main panel

# Displays the dataset


st.markdown("""
         ** :point_left: On the left you can find the parameters to set**
         """)

st.subheader(':point_down: Click RUN to generate decoys')

   
 
#%%

def decoy_fase1(loaded_smiles):
    import gc # 
    my_molecules = loaded_smiles[0].tolist()
    
    my_bar = st.progress(0)
    tamanio_total = len(my_molecules)
    t = st.empty()

    i = 0   
    
    conteo_final = []
    analisis = []
    smiles_seleccionados = []
    fp_activos = []
    fp_seleccionados = []
    smiles_x_activo_final=[]
    fp_seleccionados_x_activo_final = []
    
    for molecules in my_molecules:
        # Update progress bar.
        t.markdown("Progress: " + str(i+1) +"/" + str(tamanio_total))
        
        my_bar.progress(i + 1)
        
        conteo_1 = []               # molecules passing physicochemical filters
        filtro_tanimoto=[]          # molecules passing similarity filter        
        filtro_MCS=[]               # molecules passing MCS filter
        filtro_fw=[]                # molecules passing framework filter
        smiles_by_active= []
        fp_seleccionados_por_activo = []
        molecula_ok=molecules.strip()
        try:
            estandarizada = standardize_smiles(molecula_ok)
        except:
            st.markdown("**Oh no! There is a problem with standarization of one smiles.** :confused:")
            st.markdown("**Please check your molecule: **" + str(i+1))
            st.markdown("**That is the smiles: **" + str(molecula_ok))
            st.stop()
        i = i + 1
        nombre= "Query_" + str(i)
        my_smiles = Chem.MolFromSmiles(estandarizada) # convierte las moleculas a smiles
        core = MurckoScaffold.GetScaffoldForMol(my_smiles) # determino el framework
        framework = Chem.MolToSmiles(core) # guardo el framework en smiles
        
        #DESCRIPTORES
        MolWt = Descriptors.MolWt(my_smiles) # calculo el peso molecular
        MaxWt_DB = MolWt + lim_MW
        MinWt_DB = MolWt - lim_MW
        
        MolLogP = Descriptors.MolLogP(my_smiles) # calculo el logP
        MaxMolLogP_DB = MolLogP + lim_logP
        MinMolLogP_DB = MolLogP - lim_logP
    
        NumRotatableBonds = Descriptors.NumRotatableBonds(my_smiles) # calculo numero de enlaces rotables
        MaxNumRotatableBonds_DB = NumRotatableBonds + lim_rot_bonds
        MinNumRotatableBonds_DB = NumRotatableBonds - lim_rot_bonds
    
        NumHAcceptors = Descriptors.NumHAcceptors(my_smiles) # calculo numero de aceptores de H
        MaxNumHAcceptors_DB = NumHAcceptors + lim_N_Haccep
        MinNumHAcceptors_DB = NumHAcceptors - lim_N_Haccep
    
        NumHDonors = Descriptors.NumHDonors(my_smiles) # calculo numero de donores de H
        MaxNumHDonors_DB = NumHDonors + lim_N_Hdon
        MinNumHDonors_DB = NumHDonors - lim_N_Hdon
    
        tamanio_molec = Descriptors.HeavyAtomCount(my_smiles)
        
        # FINGERPRINT
        fp_1 = AllChem.GetMorganFingerprintAsBitVect(my_smiles,fingerprint_radio,nBits = fingerprint_lenght,useFeatures=False)
        fp_activos.append(fp_1)
            
        #BASE DE DATOS
        base_de_datos_total = os.listdir("databases")
        for base_de_datos in base_de_datos_total:
            base_de_datos = open("databases/" + base_de_datos, "r") #abro la base de datos para comparar
            for linea_DB in base_de_datos:      # le digo que repita las siguientes acciones para cada molecula en la base de datos
                linea1_DB = linea_DB.strip()
                linea2_DB = linea1_DB.split("\t")
                name=linea2_DB[0]
                print(name)
                smiles_DB = linea2_DB[1]
                framework_DB = linea2_DB[2]
                MolWt_DB = float(linea2_DB[3])
                MolLogP_DB = float(linea2_DB[4])
                NumRotatableBonds_DB = float(linea2_DB[5])
                NumHAcceptors_DB = float(linea2_DB[6])
                NumHDonors_DB = float(linea2_DB[7])
                # Physicochemical filters
                if MinWt_DB <= MolWt_DB <= MaxWt_DB:
                    if MinMolLogP_DB <= MolLogP_DB <= MaxMolLogP_DB:
                        if MinNumRotatableBonds_DB <= NumRotatableBonds_DB <= MaxNumRotatableBonds_DB:
                            if MinNumHAcceptors_DB <= NumHAcceptors_DB <= MaxNumHAcceptors_DB:
                                if MinNumHDonors_DB <= NumHDonors_DB <= MaxNumHDonors_DB:
                                    conteo_final.append(name)
                                    conteo_1.append(name)
                                    # Dissimilarity conditions                        
                                    molec_x = Chem.MolFromSmiles(smiles_DB)
                                    fp_2 = AllChem.GetMorganFingerprintAsBitVect(molec_x,fingerprint_radio,nBits = fingerprint_lenght,useFeatures=False)
                                    similarity_metric_ok = getattr(DataStructs, similarity_metric)
                                    tan_sim= similarity_metric_ok(fp_1, fp_2)
                                    if tan_sim <= float(max_similarity_limit):
                                        mols = [my_smiles,molec_x]
                                        filtro_tanimoto.append(name)
                                        res = rdFMCS.FindMCS(mols)
                                        tamanio_MCS = res.numAtoms
                                        if tamanio_MCS/tamanio_molec < lim_fraction_MCS:
                                            filtro_MCS.append(name)
                                            if framework_filter == True:
                                                if framework != framework_DB:
                                                    if not smiles_DB in smiles_seleccionados:
                                                        smiles_by_active.append(smiles_DB)
                                                        fp_seleccionados_por_activo.append(fp_2)
                                                    filtro_fw.append(name)
                                                    fp_seleccionados.append(fp_2)
                                                    smiles_seleccionados.append(smiles_DB)
                                                else:
                                                    pass
                                            if framework_filter == False:
                                                if not smiles_DB in smiles_seleccionados:
                                                    smiles_by_active.append(smiles_DB)
                                                    fp_seleccionados_por_activo.append(fp_2)
                                                filtro_fw.append(name)
                                                fp_seleccionados.append(fp_2)
                                                smiles_seleccionados.append(smiles_DB)
                                        
        base_de_datos.close()
        OK = [str(nombre), str(len(conteo_1)),str(len(filtro_tanimoto)),str(len(filtro_MCS)),str(len(filtro_fw))]
        analisis.append(OK)
        smiles_x_activo_final.append(smiles_by_active)
        fp_seleccionados_x_activo_final.append(fp_seleccionados_por_activo)
        del base_de_datos
        gc.collect()
    OK1 = pd.DataFrame(analisis)
    OK2= OK1.rename(columns={0:"Query",1:"selected by physicochemical properties",2:"pass the Tc filter",3:"pass the fMCS filter",4:"Obtained decoys"})
    resultado_fase_1 = [conteo_final, smiles_x_activo_final,fp_activos, fp_seleccionados_x_activo_final]
    
    st.success('Decoys have been successfully obtained for each loaded molecule!!! ')

    # st.write("---------------------------------------------------------------")
    st.write("Molecules that passed the physicochemical properties filters: " + str(len(conteo_final)))
    st.write("Molecules that passed the dissimilarity (structural) filters: " + str(len(smiles_seleccionados)))
    st.write("Of which: " + str(len(set(smiles_seleccionados))) + " are different")
    # st.write("if there is enough number of decoys by loaded molecule, " + str(max_decoys) +" have been randombly selected by each loaded molecule")

    return (resultado_fase_1,OK2)



#%% To take XXX random decoys from each list

def duplicates_filter(lista_resultados):
    
    fp_seleccionados_x_activo_final = lista_resultados[3]
    index_fp_no_duplicados = []
    smiles_seleccionados_ok = []
    i = 0
    random.seed(i)
    for decoys_list in lista_resultados[1]:
        fps_ok = fp_seleccionados_x_activo_final[i]
        if len(decoys_list) >= max_decoys:
            # selected_decoys = random.sample(decoys_list, max_decoys)
            selected_decoys = random.sample(list(enumerate(decoys_list)), max_decoys)
            for x in selected_decoys:
                smiles_seleccionados_ok.append(x[1])
                index_fp_no_duplicados.append(fps_ok[x[0]])
        else:
            for y,x in enumerate(decoys_list):
                smiles_seleccionados_ok.append(x)
                index_fp_no_duplicados.append(fps_ok[y])
        i = i +1 
 

#%%
    # All comparison by tanimoto
       
    comparando_todos = []
    i = 0
    
    for fp_decoy in index_fp_no_duplicados:
        smile_decoy = smiles_seleccionados_ok[i]
        coef_tan= []
        for fp_activo in lista_resultados[2]:
            tan_sim=DataStructs.TanimotoSimilarity(fp_decoy, fp_activo)
            coef_tan.append(tan_sim)
        if max(coef_tan) < max_similarity_limit_all:
            if not smile_decoy in comparando_todos:
                comparando_todos.append(smile_decoy)
        i = i + 1
    st.write("Finally, " + str(len(comparando_todos)) + " passed the Tanimoto filter by comparing all loaded molecules")
    st.success('Congratulations, you have obtained ' + str(len(comparando_todos)) + " decoys!!!")

    # st.write("---------------------------------------------------------------")
    df = pd.DataFrame(comparando_todos)
    return df


#%%
# Funcion para exportar el archivo

def filedownload(df):
    csv = df.to_csv(index=False,header=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="generated_decoys.csv">Download CSV File with your decoys</a>'
    return href

def filedownload1(df):
    csv = df.to_csv(index=False,header=True)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="decoys_analysis.csv">Download CSV File with the table</a>'
    return href

def filedownload2(df):
    csv = df.to_csv(index=False,header=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="decoys_settings.csv">Download CSV File with your settings</a>'
    return href


#%% Settings
def setting_info():
    from datetime import date
    today = date.today()
    fecha = today.strftime("%d/%m/%Y")
    settings = []
    settings.append(["Decoys generated at: " , fecha])
    settings.append(["Physicochemical features limits:",""])
    settings.append(["MW" , "+/- " +  str(lim_MW)])
    settings.append(["logP" , "+/- " + str(lim_logP)])
    settings.append(["Num_rotable_bonds" , "+/- " + str(lim_rot_bonds)])
    settings.append(["Num_H_acceptors" , "+/- " + str(lim_N_Haccep)])
    settings.append(["Num_H_donors" , "+/- " + str(lim_N_Hdon)])

    settings.append(["Topological features ---> Dissimilarty conditions",""])
    settings.append(["Morgan Fingerprints",""])
    settings.append(["Fingerprint radio:" , str(fingerprint_radio)])
    settings.append(["Fingerprint lenght:" , str(fingerprint_lenght)])
    if lim_fraction_MCS == 1.0:
        settings.append(["Not set a limit for the fraction of the Maximum Common Substructure",""])
    else:
        settings.append(["Limit of fraction of the Maximum Common Substructure: " , str(lim_fraction_MCS)])
    settings.append(["Decoys with different framework:" , str(framework_filter)])
    settings.append(["Max Tc similarity between decoys and other actives:", str(max_similarity_limit)])
    settings.append(["Max number of decoys by loaded molecule:", str(max_decoys)])
    settings_df = pd.DataFrame(settings)
    return settings_df


#%%
# ---------------------------------#

if uploaded_file_1 is not None:
    run = st.button("RUN")
    if run == True:
        loaded_smiles = pd.read_csv(uploaded_file_1,sep="\t",header=None)
        lista_resultados = decoy_fase1(loaded_smiles)
        df = duplicates_filter(lista_resultados[0])
        st.markdown(":point_down: **Here you can dowload the generated decoys**", unsafe_allow_html=True)
        st.markdown(filedownload(df), unsafe_allow_html=True)
        st.markdown(text)
        st.markdown(":point_down: **Here you can see a little analysis of the process**", unsafe_allow_html=True)
        st.write(lista_resultados[1])
        
        st.markdown(":point_down: **Here you can dowload this table in a csv file**", unsafe_allow_html=True)
        st.markdown(filedownload1(lista_resultados[1]), unsafe_allow_html=True)
        st.markdown(text)
        settings_df = setting_info()
        st.markdown(":point_down: **Here you can download your settings**", unsafe_allow_html=True)
        st.markdown(filedownload2(settings_df), unsafe_allow_html=True)

        st.balloons()

else:
    st.info('Awaiting for TXT file to be uploaded.')
    if st.button('Press to use Example Dataset'):
    #     st.markdown('The **Diabetes** dataset is used as the example.')
        st.write("We sorry, we don't have a example file yet ;)")
    #     heatmap(df_ok)






