import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import streamlit as st

# Load the CSV file into a pandas dataframe
st.set_page_config(page_title="Genomics Master")
df = pd.read_csv('clinvar_conflict.csv')
df["CLNSIGINCL"] = df["CLNSIGINCL"].replace({"B": 4, "LB": 3, "U": 2, "LP": 1, "P": 0})
df.drop(['id','Unnamed: 32'], axis=1,inplace=True)
st.title("Cancer Detection Using Classification")
st.markdown(
         f"""
         <style>
         .stApp {{
             background-image: url("https://4kwallpapers.com/images/wallpapers/macos-monterey-stock-green-dark-mode-layers-5k-6016x3384-5890.jpg");
             background-attachment: fixed;
             background-size: cover
         }}
         </style>
         """,
         unsafe_allow_html=True
     )
# Accept user inputs for the feature columns
CHROM = st.number_input('CHROM')
POS = st.number_input('POS')
REF = st.number_input('REF')
ALT = st.number_input('ALT')
AF_ESPn = st.number_input('AF_ESP')
AF_EXAC = st.number_input('AF_EXAC')
AF_TGP = st.number_input('AF_TGP')
CLNDISDB = st.number_input('CLNDISDB')
CLNDISDBINCL = st.number_input('CLNDISDBINCL')
CLNDN = st.number_input('CLNDN')
CLNDNINCL = st.number_input('CLNDNINCL')
CLNHGVS = st.number_input('CLNHGVS')
# CLNSIGINCL = st.number_input('CLNSIGINCL')
CLNVC = st.number_input('CLNVC')
CLNVI = st.number_input('CLNVI')
MC = st.number_input('MC')
ORIGIN = st.number_input('ORIGIN')
SSR = st.number_input('SSR')
CLASS = st.number_input('CLASS')
Allele = st.number_input('Allele')
Consequence = st.number_input('Consequence')
IMPACT = st.number_input('IMPACT')
SYMBOL = st.number_input('SYMBOL')
Feature_type = st.number_input('Feature_type')
Feature = st.number_input('Feature')
BIOTYPE = st.number_input('BIOTYPE')
EXON = st.number_input('EXON')
INTRON = st.number_input('INTRON')
cDNA_position = st.number_input('cDNA_position')
CDS_position = st.number_input('CDS_position')
Protein_position = st.number_input('Protein_position')
Amino_acids = st.number_input('Amino_acids')
Codons = st.number_input('Codons')
DISTANCE = st.number_input('DISTANCE')
STRAND = st.number_input('STRAND')
BAM_EDIT = st.number_input('BAM_EDIT')
SIFT = st.number_input('SIFT')
PolyPhen = st.number_input('PolyPhen')
MOTIF_NAME = st.number_input('MOTIF_NAME')
MOTIF_POS+AP1 = st.number_input('MOTIF_POS+AP1')
HIGH_INF_POS = st.number_input('HIGH_INF_POS')
MOTIF_SCORE_CHANGE = st.number_input('MOTIF_SCORE_CHANGE')
LoFtool = st.number_input('LoFtool')
CADD_PHRED = st.number_input('CADD_PHRED')
CADD_RAW+AT1 = st.number_input('CADD_RAW+AT1')
BLOSUM62 = st.number_input('BLOSUM62')



# Use the model to make predictions on the user inputs
model =LogisticRegression(C=0.1)
model.fit(df.drop('diagnosis', axis=1), df['diagnosis'])
prediction = model.predict([[radius_mean, texture_mean, perimeter_mean, area_mean,smoothness_mean, compactness_mean, concavity_mean, concave_points_mean,symmetry_mean,fractal_dimension_mean,radius_se,texture_se,perimeter_se,area_se,smoothness_se,compactness_se,concavity_se,concave_points_se,symmetry_se,fractal_dimension_se,radius_worst,texture_worst,perimeter_worst,area_worst,smoothness_worst,compactness_worst,concavity_worst,concave_points_worst,symmetry_worst,fractal_dimension_worst]])

# Display the prediction to the user
if prediction[4] == 4:
    st.write('Prediction: Benign')
if prediction[3] == 3:
     st.write('prediction: Likely Benign')
if prediction[2] == 2:
     st.write('Prediction: Uncertain Significance')
if prediction[1] == 1:
     st.write('Prediction: Likely Pathogenic')
else:
    st.write('Prediction: Pathogenic')
