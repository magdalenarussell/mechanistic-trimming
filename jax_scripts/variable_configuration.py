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

    def R_input_domain_data_path(self):
        path = self.root_path + '/meta_data/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound)
        file_name = path + '/frame_data.tsv'
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

    def validation_predictions_data_path(self, l2, validation_annotation):
        path = self.root_path + '/' + self.annotation_type + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type + '/validation_prediction'
        file_name = path + '/' + str(validation_annotation) + '_predicted_dist_data_L2' + str(l2) + '.tsv'
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
        trims = self.choice_colname.replace('_ligation-mh', '')
        trims = trims.split('_')[0]
        trims = trims.split('-')
        trims = [item + '_trim' for item in trims]
        if 'ligation-mh' in self.choice_colname:
            trims.append('ligation_mh')
        self.choice_colname = trims
        return(self.choice_colname)

    def process_group_colnames(self):
        self.group_colname = self.group_colname.split('_')[0]
        self.group_colname = self.group_colname.split('-')
        self.group_colname = [item + '_gene_group' for item in self.group_colname]
        return(self.group_colname)

    def get_all_mh_variables(self, overlap_list=[0,1,2,3,4], prop=True, pos=['up', 'mid', 'down']):
        if prop is True:
            if 'mh' in self.variable_colnames:
                self.variable_colnames.remove('mh')
            elif 'interior_mh' in self.variable_colnames:
                self.variable_colnames.remove('interior_mh')
        elif prop is False:
            if 'mh_count' in self.variable_colnames:
                self.variable_colnames.remove('mh_count')
            elif 'interior_mh_count' in self.variable_colnames:
                self.variable_colnames.remove('interior_mh_count')

        mh_vars = []
        for o in overlap_list:
            for p in pos:
                if o == 0 and p == 'mid':
                    continue
                if prop is True:
                    var = f'mh_prop_{p}_overlap_{o}'
                else:
                    var = f'mh_count_{p}_overlap_{o}'
                mh_vars.append(var)
        self.variable_colnames = self.variable_colnames + mh_vars
        return(self.variable_colnames)

    def get_all_mh_length_interaction_variables(self, overlap_list=[0,1,2,3,4], prop=True):
        assert 'mh_length_interaction' in self.variable_colnames or 'mh_count_length_interaction' in self.variable_colnames, "mh_length_interaction is not a variable column"
        if prop is True:
            self.variable_colnames.remove('mh_length_interaction')
        elif prop is False:
            self.variable_colnames.remove('mh_count_length_interaction')

        pos = ['up', 'down']
        lengths = ['j_length', 'v_length']
        mh_vars = []
        for o in overlap_list:
            for i in range(len(pos)):
                if prop is True:
                    var = f'mh_prop_{pos[i]}_overlap_{o}_{lengths[i]}_interaction'
                else:
                    var = f'mh_count_{pos[i]}_overlap_{o}_{lengths[i]}_interaction'
                mh_vars.append(var)

        self.variable_colnames = self.variable_colnames + mh_vars
        return(self.variable_colnames)

    def get_all_base_count_length_interaction_variables(self, side):
        assert side + '_base_count_length_interaction' in self.variable_colnames, "base_count_length_interaction is not a variable column"
        self.variable_colnames.remove(side + '_base_count_length_interaction')

        bases = ['AT', 'GC']
        trims = [t for t in self.choice_colname if 'trim' in t]
        base_vars = [trim + '_' + side + '_base_count_' + base + '_prop' for trim in trims for base in bases if base + side != 'AT' + side]

        interactions = []
        for var in base_vars:
            if 'v_trim' in var:
                var = var + '_v_length_interaction'
                interactions.append(var)
            elif 'j_trim' in var:
                var = var + '_j_length_interaction'
                interactions.append(var)
        self.variable_colnames = self.variable_colnames + interactions
        return(self.variable_colnames)

    def get_all_base_variables(self, side):
        assert side + '_base_count' in self.variable_colnames, "base_count is not a variable colname"
        self.variable_colnames.remove(side + '_base_count')
        bases = ['AT', 'GC']
        trims = [t for t in self.choice_colname if 'trim' in t]
        base_vars = [trim + '_' + side + '_base_count_' + base for trim in trims for base in bases if base + side != 'AT5end']
        self.variable_colnames = self.variable_colnames + base_vars
        return(self.variable_colnames)

    def get_all_base_prop_variables(self, side):
        assert side + '_base_count_prop' in self.variable_colnames, "base_count_prop is not a variable colname"
        self.variable_colnames.remove(side + '_base_count_prop')
        bases = ['AT_prop', 'GC_prop']
        trims = [t for t in self.choice_colname if 'trim' in t]
        base_vars = [trim + '_' + side + '_base_count_' + base for trim in trims for base in bases if base + side != 'AT_prop5end']
        self.variable_colnames = self.variable_colnames + base_vars
        return(self.variable_colnames)

    def get_all_motif_variables(self):
        assert 'motif' in self.variable_colnames, "motif is not a variable colname"
        assert self.left_nuc_motif_count > 0, "left motif size must be greater than zero"
        assert self.right_nuc_motif_count > 0, "right motif size must be greater than zero"
        self.variable_colnames.remove('motif')
        trims = [t for t in self.choice_colname if 'trim' in t]
        motif_vars_5 = [trim + '_motif_5end_pos' + str(pos) for trim in trims for pos in range(self.left_nuc_motif_count, 0, -1)]
        motif_vars_3 = [trim + '_motif_3end_pos' + str(pos) for trim in trims for pos in range(1, self.right_nuc_motif_count+1)]
        self.variable_colnames = self.variable_colnames + motif_vars_5 + motif_vars_3
        return(self.variable_colnames)

    def get_ligation_mh_variables(self):
        assert 'ligation_mh' in self.variable_colnames, "ligation_mh is not a varaible colname"
        return(self.variable_colnames)

    def get_all_length_variables(self):
        assert 'length' in self.variable_colnames, "length is not a variable colname"
        self.variable_colnames.remove('length')
        trims = []
        choices = [t for t in self.choice_colname if 'trim' in t]
        for c in choices:
            trims.append(c[0] + '_length')
        self.variable_colnames = self.variable_colnames + trims
        return(self.variable_colnames)

    def process_model_parameters(self):
        self.choice_colname = self.process_choice_colnames()
        self.group_colname = self.process_group_colnames()
        if 'mh' in self.variable_colnames:
            self.variable_colnames = self.get_all_mh_variables()
        if 'interior_mh' in self.variable_colnames:
            self.variable_colnames = self.get_all_mh_variables(overlap_list=[1, 2, 3, 4], pos = ['mid'])
        if 'mh_count' in self.variable_colnames:
            self.variable_colnames = self.get_all_mh_variables(prop=False)
        if 'interior_mh_count' in self.variable_colnames:
            self.variable_colnames = self.get_all_mh_variables(overlap_list=[1, 2, 3, 4], prop=False, pos = ['mid'])
        if '5end_base_count' in self.variable_colnames:
            self.variable_colnames = self.get_all_base_variables('5end')

        if '3end_base_count' in self.variable_colnames:
            self.variable_colnames = self.get_all_base_variables('3end')

        if '5end_base_count_prop' in self.variable_colnames:
            self.variable_colnames = self.get_all_base_prop_variables('5end')

        if '3end_base_count_prop' in self.variable_colnames:
            self.variable_colnames = self.get_all_base_prop_variables('3end')

        if 'motif' in self.variable_colnames:
            self.variable_colnames = self.get_all_motif_variables()

        if 'length' in self.variable_colnames:
            self.variable_colnames = self.get_all_length_variables()

        if 'mh_length_interaction' in self.variable_colnames:
            self.variable_colnames = self.get_all_mh_length_interaction_variables()
        if 'mh_count_length_interaction' in self.variable_colnames:
            self.variable_colnames = self.get_all_mh_length_interaction_variables(prop=False)
        if '3end_base_count_length_interaction' in self.variable_colnames:
            self.variable_colnames = self.get_all_base_count_length_interaction_variables('3end')
        if 'ligation_mh' in self.variable_colnames:
            self.variable_colnames = self.get_ligation_mh_variables()

        return(self)
