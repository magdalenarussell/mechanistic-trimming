class global_paramaters():
    def __init__(self, param_config, root_path, annotation_type, param_group, left_motif_size, right_motif_size, model_type):
        self.param_config = param_config
        self.trim_type = None
        self.productivity = None
        self.motif_type = None
        self.gene_name = None
        self.upper_trim_bound = None
        self.lower_trim_bound = None
        self.insertions = None
        self.model_group = None
        self.gene_weight_type = None
        self.root_path = root_path
        self.annotation_type = annotation_type
        self.param_group = param_group
        self.left_nuc_motif_count = left_motif_size
        self.right_nuc_motif_count = right_motif_size
        self.model_type = model_type

    def set(self):
        self.trim_type = getattr(self.param_config, "TRIM_TYPE")
        self.productivity = getattr(self.param_config, "PRODUCTIVITY")
        self.motif_type = getattr(self.param_config, "MOTIF_TYPE")
        self.gene_name = getattr(self.param_config, "GENE_NAME")
        self.upper_trim_bound = getattr(self.param_config, "UPPER_TRIM_BOUND")
        self.lower_trim_bound = getattr(self.param_config, "LOWER_TRIM_BOUND")
        self.insertions = getattr(self.param_config, "INSERTIONS")
        self.model_group = getattr(self.param_config, "MODEL_GROUP")
        self.gene_weight_type = getattr(self.param_config, "GENE_WEIGHT_TYPE")
        return(self)

    def R_processed_data_path(self):
        output_location = self.root_path + '/' + self.annotation_type + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/processed_data/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type + '/processed_data.tsv'
        return(output_location)

