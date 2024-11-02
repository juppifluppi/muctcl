from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdFingerprintGenerator, MACCSkeys, Descriptors, Draw
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.DataStructs import cDataStructs
from io import StringIO
from mordred import Calculator, descriptors
import numpy as np
import pandas as pd
import seaborn as sns
import sys, os, shutil
import matplotlib.pyplot as plt
import streamlit as st
from streamlit_ketcher import st_ketcher
import time
import subprocess
from PIL import Image
import uuid
from filelock import Timeout, FileLock
from streamlit_echarts import st_echarts
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from scopy.ScoPretreat import pretreat
import scopy.ScoDruglikeness
from dimorphite_dl import DimorphiteDL
import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

def cooling_highlight(val):
   color = "red" if val < 50 else "green"                    
   return f'background-color: {color}'

calc = Calculator(descriptors, ignore_3D=False)

def fingerprint_rdkit(self,rad,fplength) -> np.ndarray:
    fp_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=rad,fpSize=fplength)
    return fp_gen.GetCountFingerprintAsNumPy(self).astype(int)
   
def standardize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    clean_mol = rdMolStandardize.Cleanup(mol) 
    parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol)
    uncharger = rdMolStandardize.Uncharger()
    uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)    
    te = rdMolStandardize.TautomerEnumerator()
    taut_uncharged_parent_clean_mol = te.Canonicalize(uncharged_parent_clean_mol)     
    taut_uncharged_parent_clean_mol_addH = Chem.AddHs(taut_uncharged_parent_clean_mol)
    Chem.SanitizeMol(taut_uncharged_parent_clean_mol_addH)
    return taut_uncharged_parent_clean_mol_addH

def cooling_highlight(val):
   color = 'green' if val > 49 else "red"                    
   return f'background-color: {color}'

st.title('BILE/MUCIN interaction model')

with st.form(key='my_form_to_submit'):
    with st.expander("More information"):
        
        st.caption(""":black[Background]""")
        st.caption("""muctcl predicts interactions of drugs with bile and mucin, in binary and ternary mixtures. Models are trained based on classifications derived from ¹H-NMR measurements. Binomial logistic regression models are used for classifications
        of binary drug-bile and drug-mucin mixtures, whereas a 4-class random forest model is used for a prediction of mixtures.""")
        
        st.caption("""The software is hosted at our [github page](https://github.com/juppifluppi/muctcl).""")
 
        st.caption("""Version 1.0 (02.11.2024)""")
 
    SMI = st.text_input('Enter [SMILES code](https://pubchem.ncbi.nlm.nih.gov//edit3/index.html) of drug to load', '') 
    
#    on = st.toggle('Use drawn structure')
#    with st.expander("SMILES editor"):
#        drawer = st_ketcher()
#        st.caption("Click on Apply to save the drawn structure as input.")
#
#    if on:
#        SMI=drawer


    on3 = st.toggle('Perform batch calculation',key="16")    
    with st.expander("Batch settings"):
        col1, col2 = st.columns(2)
        with col1:
            st.write("Names of compounds")
            NAMESx = st.text_area(label="Input names of compounds separated by linebreaks",key="17")
        with col2:
            st.write("SMILES codes")
            SMILESx = st.text_area(label="Input SMILES of compounds separated by linebreaks",key="18")

    
    emoji = ''
    label = ' Predict'    
    submit_button = st.form_submit_button(label=f'{emoji} {label}')

if submit_button:
# try:
    lock = FileLock("lockfile.lock")
    with st.spinner('WAITING IN QUEUE ...'):
        try:
            lock.acquire(timeout=20)
        except Timeout:
            os.remove("lockfile.lock")
            lock = FileLock("lockfile.lock")
            lock.acquire()

    try:
        if on3 is False:
    
            for es in ["descriptors.csv","results.csv","results3.csv","results4.csv","results5.csv","results2.csv"]:
                try:
                    os.remove(es)
                except:
                    pass
    
            SMI=SMI
                                 
            mol = standardize(SMI)
            maccskeys = MACCSkeys.GenMACCSKeys(mol)     
            rdk5fp1 = fingerprint_rdkit(mol,5,2048)

            try:
                sdm = pretreat.StandardizeMol()
                molx = sdm.disconnect_metals(mol)    
                logd = scopy.ScoDruglikeness.molproperty.CalculateLogD(molx)
                mr = scopy.ScoDruglikeness.molproperty.CalculateMolMR(molx)    
                tcl1 = ( ( logd - 1.510648) / 1.708574 ) * 1.706694
                tcl2 = ( ( mr - 90.62889 ) / 35.36033 ) * 2.4925333    
                tcl3 = 1 / ( 1 + ( 2.718281828459045 ** ( -1 * ( 0.9872289 + tcl1 + tcl2 ) ) ) )    
                #st.write("logD: " + str(round(logd,2)))
                #st.write("CrippenMR: " + str(round(mr,2)))
                #st.write("TC/L interaction probability: " + str(int(round(tcl3*100,2))) + " %")
                #st_echarts(liquidfill_option)

            except:
                st.write("Something is wrong with your SMILES code.")
                st.stop()
           
            with open("descriptors.csv","a") as f:
                for o in range(0,len(maccskeys)):
                    f.write("maccs_"+str(o)+"\t")    
                for o in range(0,len(rdk5fp1)):
                    f.write("rdk5fp1_"+str(o)+"\t")      
                f.write("\n")  
                
            with open("descriptors.csv","a") as f:
                for o in range(0,len(maccskeys)):
                    f.write(str(maccskeys[o])+"\t") 
                for o in range(0,len(rdk5fp1)):
                    f.write(str(rdk5fp1[o])+"\t")   
                f.write("\n")
                                                                                     
            process3=subprocess.Popen(["Rscript", "predict.R"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            process3.communicate()
           
            df2 = pd.read_csv(r'results.csv')
            df3 = pd.read_csv(r'results2.csv')    
            df4 = pd.read_csv(r'results3.csv') 
            df5 = pd.read_csv(r'results4.csv') 
            df6 = pd.read_csv(r'results5.csv') 
           
            col1, col2 = st.columns(2)

            with col1: 
                style_heading = f"font-size: 15px; text-align: center;"
                st.markdown(f"<h1 style='{style_heading}'>MOLECULE</h1>", unsafe_allow_html=True)
                im = Draw.MolToImage(Chem.MolFromSmiles(SMI),fitImage=True)
                st.image(im)
                st.markdown(f"<h1 style='{style_heading}'>BILE PREDICTION</h1>", unsafe_allow_html=True)
                liquidfill_option = {
                "series": [{"type": "liquidFill", "data": [int(round(tcl3*100,2))]}]
                }
                st_echarts(liquidfill_option,key="23456")
                #st.image("tc.png",width=80)
                st.markdown("<img src='https://raw.githubusercontent.com/juppifluppi/muctcl/refs/heads/main/tc.png' width='150' style='display: block; margin: 0 auto;'>" , unsafe_allow_html=True)
               
            with col2: 
                st.markdown(f"<h1 style='{style_heading}'>MIXTURE PREDICTION</h1>", unsafe_allow_html=True)
                options = {
                    "tooltip": {"trigger": "item"},
                    "legend": {"top": "5%", "left": "center"},
                    "series": [
                        {
                            "name": "probability [%]",
                            "type": "pie",
                            "radius": ["10%", "30%"],
                            "avoidLabelOverlap": False,
                            "itemStyle": {
                                "borderRadius": 2,
                                "borderColor": "#fff",
                                "borderWidth": 2,
                            },
                            "label": {"show": False, "position": "center"},
                            "emphasis": {
                                "label": {"show": True, "fontSize": "20", "fontWeight": "bold"}
                            },
                            "labelLine": {"show": False},
                            "data": [
                                {"value": int(df3.iloc[0, 0]*100), "name": "bile+mucin interacting", "itemStyle": {"color": "#ff9999"}},
                                {"value": int(df4.iloc[0, 0]*100), "name": "mucin interacting", "itemStyle": {"color": "#99ff99"}},
                                {"value": int(df5.iloc[0, 0]*100), "name": "bile interacting", "itemStyle": {"color": "#66b3ff"}},
                                {"value": int(df6.iloc[0, 0]*100), "name": "non-interacting", "itemStyle": {"color": "#ffcc99"}},
                            ],
                        }
                    ],
                }
                st_echarts(
                    options=options, height="290px",
                )
                st.markdown(f"<h1 style='{style_heading}'>MUCIN PREDICTION</h1>", unsafe_allow_html=True)
                liquidfill_option2 = {
                "series": [{"type": "liquidFill", "data": [int(df2.iloc[0, 0]*1)]}]
                }
                st_echarts(liquidfill_option2,key="3456")                                                         
                #st.image("muc2.png",width=80)
                st.markdown("<img src='https://raw.githubusercontent.com/juppifluppi/muctcl/refs/heads/main/muc2.png' width='150' style='display: block; margin: 0 auto;'>" , unsafe_allow_html=True)

            for es in ["descriptors.csv","results.csv","results2.csv","results3.csv","results4.csv","results5.csv"]:
                try:
                    os.remove(es)
                except:
                    pass
   
        if on3 is True:

           SMILESx=SMILESx.split('\n')           
           NAMESx=NAMESx.split('\n')  

           for es in ["descriptors.csv","results.csv","results3.csv","results4.csv","results5.csv","results2.csv"]:
                try:
                    os.remove(es)
                except:
                    pass
           omoo=[]

           for cx in range(0,len(SMILESx)):
               SMI=SMILESx[cx]
                          
               mol = standardize(SMI)
               maccskeys = MACCSkeys.GenMACCSKeys(mol)     
               rdk5fp1 = fingerprint_rdkit(mol,5,2048)

               try:
                   sdm = pretreat.StandardizeMol()
                   molx = sdm.disconnect_metals(mol)    
                   logd = scopy.ScoDruglikeness.molproperty.CalculateLogD(molx)
                   mr = scopy.ScoDruglikeness.molproperty.CalculateMolMR(molx)    
                   tcl1 = ( ( logd - 1.510648) / 1.708574 ) * 1.706694
                   tcl2 = ( ( mr - 90.62889 ) / 35.36033 ) * 2.4925333    
                   tcl3 = 1 / ( 1 + ( 2.718281828459045 ** ( -1 * ( 0.9872289 + tcl1 + tcl2 ) ) ) )     
                   omoo.append(int(round(tcl3*100,2)))
    
               if cx == 0:

                   with open("descriptors.csv","a") as f:
                       for o in range(0,len(maccskeys)):
                           f.write("maccs_"+str(o)+"\t")    
                       for o in range(0,len(rdk5fp1)):
                           f.write("rdk5fp1_"+str(o)+"\t")      
                       f.write("\n")  
                
               with open("descriptors.csv","a") as f:
                   for o in range(0,len(maccskeys)):
                       f.write(str(maccskeys[o])+"\t") 
                   for o in range(0,len(rdk5fp1)):
                       f.write(str(rdk5fp1[o])+"\t")   
                   f.write("\n")
            
           process3=subprocess.Popen(["Rscript", "predict.R"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
           process3.communicate()
                                           
           df2 = pd.read_csv(r'results.csv')
           df3 = pd.read_csv(r'results2.csv')        
           df4 = pd.read_csv(r'results3.csv') 
           df5 = pd.read_csv(r'results4.csv') 
           df6 = pd.read_csv(r'results5.csv') 
           
           dfx = pd.DataFrame(columns=['Compound', "Mixture 1", "Mixture 2", "Mixture 3", "Mixture 4", "Isolated 1", "Isolated 2"])
           dfx["Compound"]=NAMESx
           dfx["Mixture 1"]=(df3.iloc[:, 0].astype(float))*100
           dfx["Mixture 1"]=dfx.iloc[:, 1].astype(int)
           dfx["Mixture 2"]=(df4.iloc[:, 0].astype(float))*100
           dfx["Mixture 2"]=dfx.iloc[:, 1].astype(int)
           dfx["Mixture 3"]=(df5.iloc[:, 0].astype(float))*100
           dfx["Mixture 3"]=dfx.iloc[:, 1].astype(int)
           dfx["Mixture 4"]=(df6.iloc[:, 0].astype(float))*100
           dfx["Mixture 4"]=dfx.iloc[:, 1].astype(int)
           dfx["Isolated 1"]=(df2.iloc[:, 0].astype(float))*100
           dfx["Isolated 1"]=dfx.iloc[:, 1].astype(int)
           dfx["Isolated 2"]=(omoo)
           dfx["Isolated 2"]=dfx.iloc[:, 1].astype(int)
    
           #dfx.reset_index(inplace=True)               
           st.dataframe(dfx.style.applymap(cooling_highlight,subset=["Mixture 1", "Mixture 2"]))    
    
    finally:
        lock.release()

#    except:
#        st.write("Something went wrong. Cannot parse molecules! Please verify your structures.")  

