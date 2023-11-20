import os

class global_paramaters():
    def __init__(self, param_config, root_path, project_path, annotation_type, param_group, left_motif_size, right_motif_size, model_type):
        self.param_config = param_config
        self.root_path = root_path
        self.project_path = project_path
        self.annotation_type = annotation_type
        self.param_group = param_group
        self.left_nuc_motif_count = left_motif_size
        self.right_nuc_motif_count = right_motif_size
        self.model_type = model_type
        self.trim_type = getattr(self.param_config, "TRIM_TYPE")
        self.productivity = getattr(self.param_config, "PRODUCTIVITY")
        self.motif_type = getattr(self.param_config, "MOTIF_TYPE")
        self.gene_name = getattr(self.param_config, "GENE_NAME")
        self.upper_trim_bound = getattr(self.param_config, "UPPER_TRIM_BOUND")
        self.lower_trim_bound = getattr(self.param_config, "LOWER_TRIM_BOUND")
        self.insertions = getattr(self.param_config, "INSERTIONS")
        self.model_group = getattr(self.param_config, "MODEL_GROUP")
        self.gene_weight_type = getattr(self.param_config, "GENE_WEIGHT_TYPE")

    def R_processed_data_path(self, annotation = None):
        if annotation == None:
            annotation = self.annotation_type
        path = self.root_path + '/' + annotation + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type
        file_name = path + '/processed_data.tsv'
        return(file_name)

    def R_subsampling_processed_data_path(self, prop, annotation = None):
        if annotation == None:
            annotation = self.annotation_type
        path = self.root_path + '/' + annotation + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type +'/temp_subsampling_exp/prop' + prop
        return(path)

    def model_output_path(self, l2):
        path = self.project_path + '/trained_models/' + self.annotation_type + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type
        os.makedirs(path, exist_ok=True)
        file_name = path + '/trained_model_L2' + str(l2) + '.pkl'
        return(file_name)

    def predictions_data_path(self, l2):
        path = self.root_path + '/' + self.annotation_type + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type
        file_name = path + '/predicted_dist_data_L2' + str(l2) + '.tsv'
        return(file_name)

    def trained_coefs_path(self, l2):
        path = self.root_path + '/' + self.annotation_type + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type
        file_name = path + '/trained_coefs_L2' + str(l2) + '.tsv'
        return(file_name)

    def subsampling_coefs_path(self, prop, l2):
        path = self.root_path + '/' + self.annotation_type + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type + '/subsampling_experiment'
        os.makedirs(path, exist_ok=True)
        file_name = path + '/data_prop' + str(prop) + '_coefs_L2' + str(l2) + '.tsv'
        return(file_name)

    def model_eval_results_path(self, l2):
        path = self.root_path + '/' + self.annotation_type + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type
        file_name = path + '/model_evaluation_results_L2' + str(l2) + '.tsv'
        return(file_name)


class model_specific_parameters():
    def __init__(self, param_config, model_type_config, left_motif_size, right_motif_size):
        self.param_config = param_config
        self.model_type_config = model_type_config
        self.left_nuc_motif_count = left_motif_size
        self.right_nuc_motif_count = right_motif_size
        self.variable_colnames = getattr(self.model_type_config, "VARIABLE_COLNAMES")
        self.count_colname = getattr(self.model_type_config, "COUNT_COLNAME")
        self.group_colname = getattr(self.param_config, "GENE_NAME")
        self.repeat_obs_colname = getattr(self.param_config, "REPEAT_OBS_COLNAME")
        self.choice_colname = getattr(self.param_config, "TRIM_TYPE")

    def process_choice_colnames(self):
        self.choice_colname = self.choice_colname.split('_')[0]
        self.choice_colname = self.choice_colname.split('-')
        self.choice_colname = [item + '_trim' for item in self.choice_colname]
        return(self.choice_colname)

    def process_group_colnames(self):
        self.group_colname = self.group_colname.split('_')[0]
        self.group_colname = self.group_colname.split('-')
        self.group_colname = [item + '_gene_group' for item in self.group_colname]
        return(self.group_colname)

    def get_all_mh_variables(self, overlap_list=[0,1,2,3,4]):
        assert 'mh' in self.variable_colnames, "mh is not a variable colname"
        self.variable_colnames.remove('mh')
        pos = ['up', 'mid', 'down']
        mh_vars = []
        for o in overlap_list:
            for p in pos:
                if o == 0 and p == 'mid':
                    continue
                var = f'mh_prop_{p}_overlap_{o}'
                mh_vars.append(var)
        self.variable_colnames = self.variable_colnames + mh_vars
        return(self.variable_colnames)

    def get_all_base_variables(self):
        assert 'base_count' in self.variable_colnames, "base_count is not a variable colname"
        self.variable_colnames.remove('base_count')
        sides = ['5end', '3end']
        bases = ['AT', 'GC']
        trims = self.choice_colname
        base_vars = [trim + '_' + side + '_base_count_' + base for trim in trims for side in sides for base in bases if base + side != 'AT5end']
        self.variable_colnames = self.variable_colnames + base_vars
        return(self.variable_colnames)

    def get_all_base_prop_variables(self):
        assert 'base_count_prop' in self.variable_colnames, "base_count_prop is not a variable colname"
        self.variable_colnames.remove('base_count_prop')
        sides = ['5end', '3end']
        bases = ['AT_prop', 'GC_prop']
        trims = self.choice_colname
        base_vars = [trim + '_' + side + '_base_count_' + base for trim in trims for side in sides for base in bases if base + side != 'AT_prop5end']
        self.variable_colnames = self.variable_colnames + base_vars
        return(self.variable_colnames)

    def get_all_motif_variables(self):
        assert 'motif' in self.variable_colnames, "motif is not a variable colname"
        assert self.left_nuc_motif_count > 0, "left motif size must be greater than zero"
        assert self.right_nuc_motif_count > 0, "right motif size must be greater than zero"
        self.variable_colnames.remove('motif')
        trims = self.choice_colname
        motif_vars_5 = [trim + '_motif_5end_pos' + str(pos) for trim in trims for pos in range(self.left_nuc_motif_count, 0, -1)]
        motif_vars_3 = [trim + '_motif_3end_pos' + str(pos) for trim in trims for pos in range(1, self.right_nuc_motif_count+1)]
        self.variable_colnames = self.variable_colnames + motif_vars_5 + motif_vars_3
        return(self.variable_colnames)

    def get_all_length_variables(self):
        assert 'length' in self.variable_colnames, "length is not a variable colname"
        self.variable_colnames.remove('length')
        trims = self.choice_colname
        self.variable_colnames = self.variable_colnames + trims
        return(self.variable_colnames)

    def process_model_parameters(self):
        self.choice_colname = self.process_choice_colnames()
        self.group_colname = self.process_group_colnames()
        if 'mh' in self.variable_colnames:
            self.variable_colnames = self.get_all_mh_variables()
        if 'base_count' in self.variable_colnames:
            self.variable_colnames = self.get_all_base_variables()
        if 'base_count_prop' in self.variable_colnames:
            self.variable_colnames = self.get_all_base_prop_variables()
        if 'motif' in self.variable_colnames:
            self.variable_colnames = self.get_all_motif_variables()
        if 'length' in self.variable_colnames:
            self.variable_colnames = self.get_all_length_variables()
        return(self)
