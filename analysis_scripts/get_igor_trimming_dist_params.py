import pygor3 as p3
import sys
import pandas as pd

MOD_OUTPUT_PATH=sys.argv[1]

mdl_hb = p3.get_default_IgorModel("human", "tcr_alpha")

vgenes = pd.DataFrame()

for i in range(len(mdl_hb['v_3_del'])):
    item = mdl_hb['v_3_del'][i]
    temp_data = {'v_gene': item.lbl__v_choice.values,'v_trim': item.lbl__v_3_del.values,'v_trim_prob': item.values}
    temp = pd.DataFrame(temp_data)
    vgenes = pd.concat([vgenes, temp], ignore_index = True)

jgenes = pd.DataFrame()

for i in range(len(mdl_hb['j_5_del'])):
    item = mdl_hb['j_5_del'][i]
    temp_data = {'j_gene': item.lbl__j_choice.values,'j_trim': item.lbl__j_5_del.values,'j_trim_prob': item.values}
    temp = pd.DataFrame(temp_data)
    jgenes = pd.concat([jgenes, temp], ignore_index = True)

# Create a key column in both DataFrames for the merge operation
vgenes['key'] = 1
jgenes['key'] = 1

# Perform the Cartesian join
tog = pd.merge(vgenes, jgenes, on='key')

# Drop the key column as it's no longer needed
tog.drop('key', axis=1, inplace=True)

path=MOD_OUTPUT_PATH + '/meta_data/igor_alpha_trimming_params.tsv'
tog.to_csv(path, sep='\t', index=False)
print(path)
