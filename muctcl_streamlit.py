from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdFingerprintGenerator, MACCSkeys, Descriptors, Draw
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.DataStructs import cDataStructs
from io import StringIO
from mordred import Calculator, descriptors
from rdkit.Chem import rdDepictor
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
#from dimorphite_dl import DimorphiteDL
from dimorphite_dl import protonate_smiles
import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from io import BytesIO

def show_mol(d2d,mol,legend='',highlightAtoms=[]):
    d2d.DrawMolecule(mol,legend=legend, highlightAtoms=highlightAtoms)
    d2d.FinishDrawing()
    bio = BytesIO(d2d.GetDrawingText())
    return Image.open(bio)
def show_images(imgs,buffer=5):
    height = 0
    width = 0
    for img in imgs:
        height = max(height,img.height)
        width += img.width
    width += buffer*(len(imgs)-1)
    res = Image.new("RGBA",(width,height))
    x = 0
    for img in imgs:
        res.paste(img,(x,0))
        x += img.width + buffer
    return res

#def run_dimorphite_dl(SMI):
#    result = subprocess.run(
#        ["dimorphite_dl", "--ph_min", "6.4", "--ph_max", "6.6", "--precision", "0.1", "--max_variants", "1", "--silent"],
#        input=SMI.encode(),  # pass SMI as input
#        stdout=subprocess.PIPE,
#        stderr=subprocess.PIPE,
#    )
#    if result.returncode != 0:
#        raise RuntimeError(f"Dimorphite-DL failed: {result.stderr.decode()}")
#    return result.stdout.decode().strip()

#dimorphite_dl = DimorphiteDL(
#    min_ph = 6.4,
#    max_ph = 6.6,
#    max_variants = 1,
#    label_states = False,
#    pka_precision = 0.1
#)
#dimorphite_dl: list[str] = protonate_smiles(
#    "CCC(=O)O", ph_min=6.4, ph_max=6.6, precision=0.1, max_variants=1, label_states=False
#)
#print(f"Protonated 'CCC(=O)O': {dimorphite_dl}")

def run_dimorphite(SMI):
    command = [
        "dimorphite_dl",
        "--ph_min", "6.4",
        "--ph_max", "6.6",
        "--precision", "0.1",
        "--max_variants", "1",
        "--silent",
        SMI  # This is the required positional argument
    ]
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Dimorphite-DL failed:\n{result.stderr}")
    
    return result.stdout.strip()

def cooling_highlight(val):
   color = "red" if val < 50 else "green"                    
   return f'background-color: {color}'
def highlight_max(s):
    is_max = s == s.max() # Create a boolean Series where the max value is True
    return ['background-color: green' if v else '' for v in is_max]

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

st.title('bile/mucin interaction model')

with st.form(key='my_form_to_submit'):
    with st.expander("More information"):
        
        st.caption(""":black[Background]""")
        st.caption("""muctcl predicts interactions of drugs with bile and mucin, in binary and ternary mixtures. Models are trained based on classifications derived from ¹H-NMR measurements. Binomial logistic regression models are used for classifications
        of binary drug-bile and drug-mucin mixtures, whereas a 4-class random forest model is used for a prediction of mixtures. When selecting batch mode, users can input names and corresponding SMILES codes of multiple molecules in order to output all predictions
        in an exportable table format.""")
        
        st.caption("""The software is hosted at our [github page](https://github.com/juppifluppi/muctcl).""")
 
        st.caption("""Version 1.0 (02.11.2024)""")
 
    SMI = st.text_input('Enter [SMILES code](https://pubchem.ncbi.nlm.nih.gov//edit3/index.html) of drug to predict', '') 

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
                #SMIy = str(dimorphite_dl.protonate(SMI)[0])
                SMIy = run_dimorphite(SMI)
                moly = Chem.MolFromSmiles(SMIy)
                sdm = pretreat.StandardizeMol()
                moly = sdm.disconnect_metals(moly)
                logd = scopy.ScoDruglikeness.molproperty.CalculateLogD(moly)
                mr = scopy.ScoDruglikeness.molproperty.CalculateMolMR(moly)    
                tcl1 = ( ( logd - 1.510648) / 1.708574 ) * 1.706694
                tcl2 = ( ( mr - 90.62889 ) / 35.36033 ) * 2.4925333    
                tcl3 = 1 / ( 1 + ( 2.718281828459045 ** ( -1 * ( 0.9872289 + tcl1 + tcl2 ) ) ) )   

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
                imgs = []
                im = Draw.MolDraw2DCairo(350,300)
                dopts = im.drawOptions()
                dopts.setBackgroundColour((0,0,0,0))
                imgs.append(show_mol(im,(Chem.MolFromSmiles(SMI))))
                st.image(imgs[0])
                st.markdown(f"<h1 style='{style_heading}'>BILE PREDICTION</h1>", unsafe_allow_html=True)
                liquidfill_option = {
                "series": [{"type": "liquidFill", "data": [round(tcl3*10/10,2)]}]
                }
                st_echarts(liquidfill_option,key="23456")
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
                                "borderRadius": 1,
                                "borderColor": "#000",
                                "borderWidth": 1,
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
                    options=options, height="288px",
                )
                st.markdown(f"<h1 style='{style_heading}'>MUCIN PREDICTION</h1>", unsafe_allow_html=True)
                liquidfill_option2 = {
                "series": [{"type": "liquidFill", "data": [int(df2.iloc[0, 0]*100)/100]}]
                }
                st_echarts(liquidfill_option2,key="3456")                                                         
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

               #SMIy = str(dimorphite_dl.protonate(SMI)[0])
               SMIy = run_dimorphite(SMI)
               moly = Chem.MolFromSmiles(SMIy)
               sdm = pretreat.StandardizeMol()
               moly = sdm.disconnect_metals(moly)
               logd = scopy.ScoDruglikeness.molproperty.CalculateLogD(moly)
               mr = scopy.ScoDruglikeness.molproperty.CalculateMolMR(moly)    
               tcl1 = ( ( logd - 1.510648) / 1.708574 ) * 1.706694
               tcl2 = ( ( mr - 90.62889 ) / 35.36033 ) * 2.4925333    
               tcl3 = 1 / ( 1 + ( 2.718281828459045 ** ( -1 * ( 0.9872289 + tcl1 + tcl2 ) ) ) )     
               omoo.append(round(tcl3*100,2))
    
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
           
           dfx = pd.DataFrame(columns=['Compound', "I", "II", "III", "IV", "bile", "mucin"])
           dfx["Compound"]=NAMESx
           dfx["I"]=(df3.iloc[:, 0].astype(float))*100
           dfx["I"]=dfx.iloc[:, 1].astype(int)
           dfx["II"]=(df4.iloc[:, 0].astype(float))*100
           dfx["II"]=dfx.iloc[:, 2].astype(int)
           dfx["III"]=(df5.iloc[:, 0].astype(float))*100
           dfx["III"]=dfx.iloc[:, 3].astype(int)
           dfx["IV"]=(df6.iloc[:, 0].astype(float))*100
           dfx["IV"]=dfx.iloc[:, 4].astype(int)
           dfx["bile"]=(omoo)
           dfx["bile"]=dfx.iloc[:, 5].astype(int)
           dfx["mucin"]=(df2.iloc[:, 0].astype(float))*100
           dfx["mucin"]=dfx.iloc[:, 6].astype(int)
    
           #dfx.reset_index(inplace=True)               
           st.caption("The following table gives the predictions (in %) of the 4-class model for each class (I = bile+mucin interacting; II = mucin interacting; III = bile interacting; IV = non-interacting). After that, the prediction for interaction with bile and mucin for isolated measurements are given.")
           st.dataframe(dfx.style.applymap(cooling_highlight,subset=["bile", "mucin"]))    

    finally:
        lock.release()

#    except:
#        st.write("Something went wrong. Cannot parse molecules! Please verify your structures.")  

