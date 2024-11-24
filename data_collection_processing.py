import pandas as pd
from chembl_webresource_client.new_client import new_client

#target search for SARS-cov replicase

target = new_client.target
target_query = target.search('SARS coronavirus')
target = pd.DataFrame.from_dict(target_query)

#select a protein molecule

selected_target = target.target_chembl_id[0]

#retrieve bioactivity data for SARS coronavirus 3C-like proteinase as CSV 

activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type='IC50')

df = pd.DataFrame.from_dict(res)

#check for any missing NA values:

df_cleaned = df[df.standard_value.notna()]

#convert bioactivity data to a csv file

df_cleaned.to_csv('Bioactivity_data.csv', index = False)

#labelling compounds as either active, inactive or intermediate

def lab_compounds(dataframe):

    bioactivity = []

    for value in dataframe.standard_value:
        if float(value) <= 1000:
            bioactivity.append('active')
        elif float(value) >= 10000:
            bioactivity.append('inactive')
        elif float(value) > 1000 and float(value) < 10000:
            bioactivity.append('intermediate')
    
    dataframe['Bioactivity_class'] = bioactivity 

    return dataframe

lab_compounds(df_cleaned)

#creating a preprocessed dataframe:

def preprocessing(dataframe):

  mol_id = []
  canonical_smiles = []
  standard_value = []

  for index in dataframe.molecule_chembl_id:
    mol_id.append(index)

  for index in dataframe.canonical_smiles:
    canonical_smiles.append(index)

  for index in dataframe.standard_value:
    standard_value.append(index)

  return mol_id, canonical_smiles, standard_value


mol_id, canonical_smiles, standard_value = preprocessing(df_cleaned)

bioactivity_class = list(df_cleaned.Bioactivity_class)

data_dic = {"molecule_chembl_id": mol_id, "canonical_smiles": canonical_smiles, "standard_value": standard_value, "bioactivity_class": bioactivity_class}

bioactivity_preprocessing = pd.DataFrame(data_dic)

bioactivity_preprocessing.to_csv("bioactivity_preproccessing.csv", index=False)





