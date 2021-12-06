import pygor3 as p3
import sys
import pandas as pd

output_directory = sys.argv[1]
final_output_directory = sys.argv[2]

scenarios = output_directory + "/foo_output/sampled_scenarios.csv"
print(scenarios)

default_human_model = p3.get_default_IgorModel("human", "tcr_beta")

df_scenarios = default_human_model.get_dataframe_scenarios(scenarios)

# retrieve actual junction feature quantities
junction_variables = {'v_trim':'v_3_del', 'd0_trim':'d_5_del', 'd1_trim':'d_3_del', 'j_trim':'j_5_del', 'dj_insert':'dj_ins', 'vd_insert':'vd_ins'}

for var in junction_variables.keys():
    df_scenarios[var] = default_human_model.get_realization_value_from_df_scenarios(df_scenarios, junction_variables[var])

# retrieve actual gene feature values
def get_vgene(scenario):
    v_choice = default_human_model.realization(scenario, 'v_choice')
    return v_choice

def get_dgene(scenario):
    d_choice = default_human_model.realization(scenario, 'd_gene')
    return d_choice

def get_jgene(scenario):
    j_choice = default_human_model.realization(scenario, 'j_choice')
    return j_choice

df_scenarios['v_gene_call'] = default_human_model.get_observable_from_df_scenarios(get_vgene, df_scenarios = df_scenarios)

df_scenarios['d_gene_call'] = default_human_model.get_observable_from_df_scenarios(get_dgene, df_scenarios = df_scenarios)

df_scenarios['j_gene_call'] = default_human_model.get_observable_from_df_scenarios(get_jgene, df_scenarios = df_scenarios)

vd_nucs = default_human_model.get_df_realizations_dinucl(df_scenarios, 'vd_dinucl')
df_scenarios['vd_insert_nucs'] = vd_nucs.value.apply(lambda x: "".join(x))

dj_nucs = default_human_model.get_df_realizations_dinucl(df_scenarios, 'dj_dinucl')
df_scenarios['dj_insert_nucs'] = dj_nucs.value.apply(lambda x: "".join(x))

file_name = output_directory + "/foo_output/parsed_sampled_output.csv"
df_scenarios.to_csv(file_name)
