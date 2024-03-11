import pandas as pd
import numpy as np
import jax
import jax.numpy as jnp
import jaxopt
import patsy
import dill
import pickle
from pandarallel import pandarallel
from patsy.contrasts import Sum
from sklearn.model_selection import GroupKFold

class DataPreprocessor():
    """
    A class for preprocessing and transforming data.

    Args:
        training_df (pd.DataFrame): The training DataFrame.
        variable_colnames (list): List of variable column names.
        count_colname (str): The column name for the count variable.
        group_colname (str): The column name representing group identifiers.
        repeat_obs_colname (str): The column name representing repeated observation identifiers.
        choice_colname (str): The column name representing choice identifiers.

    Attributes:
        training_df (pd.DataFrame): The training DataFrame.
        variable_colnames (list): List of variable column names.
        count_colname (str): The column name for the count variable.
        group_colname (str): The column name representing group identifiers.
        repeat_obs_colname (str): The column name representing repeated observation identifiers.
        choice_colname (str): The column name representing choice identifiers.

    Methods:
        get_mapping_dict(training_df, col):
            Create a mapping dictionary for unique values in a column.

        transform_categorical_response_vars(training_df, col, new_col):
            Transform a categorical variable into integers.

        get_contrast_matrix(training_df, col):
            Get the contrast matrix for categorical variables.

        get_dropped_contrast_var(training_df, original_col):
            Get the names of dropped contrast variables after transformation.

        transform_categorical_vars(training_df, col):
            Transform categorical variable columns into contrast columns.

        expand_multivariable(training_df, col, new_col):
            Expand multi-variable columns into a single string column.
    """
    def __init__(self, training_df, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params):
        self.training_df = training_df
        self.variable_colnames = variable_colnames
        self.choice1_variable_colnames = choice1_variable_colnames
        self.choice2_variable_colnames = choice2_variable_colnames
        self.count_colname = count_colname
        self.group_colname = group_colname
        self.repeat_obs_colname = repeat_obs_colname
        self.choice_colname = choice_colname
        self.params = params

    def get_mapping_dict(self, training_df, col):
        """
        Create a mapping dictionary for unique values in a column.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            col (str): The name of the column for which the mapping is generated.

        Returns:
            dict: A dictionary mapping unique values to their corresponding indices.
        """
        # Create a mapping dictionary for unique values in a column
        unique = sorted(list(pd.unique(training_df[col])))
        mapping = {key: index for index, key in enumerate(unique)}
        return mapping

    def transform_categorical_response_vars(self, training_df, col, new_col):
        """
        Transform a categorical variable into integers.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            col (str): The name of the column to be transformed.
            new_col (str): The name of the new integer-encoded column.

        Returns:
            pd.DataFrame: The DataFrame with the integer-encoded column added.
        """
        # Transform a categorical variable into integers
        col = getattr(self, col)
        if training_df[col].iloc[1] != int:
            new_col_name = col + '_int'
            setattr(self, new_col, new_col_name)
            var_mapping = self.get_mapping_dict(training_df, col)
            training_df[new_col_name] = training_df[col].map(var_mapping)
        return training_df

    def get_contrast_matrix(self, training_df, col, pretrain=True):
        """
        Get the contrast matrix for categorical variables.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            col (str): The name of the categorical column.

        Returns:
            pd.DataFrame: The contrast matrix for the categorical column.
        """
        first = training_df[col].iloc[0]
        if not isinstance(first, int) and not isinstance(first, float) and not isinstance(first, np.int64):
            # create sum contrasts matrix
            unique = sorted(list(pd.unique(training_df[col])))
            if not pretrain:
                unique = ['-', 'A', 'C', 'G', 'T']

            contrast = Sum().code_without_intercept(unique).matrix

            # Create an identity matrix with the number of categories
            identity_matrix = np.eye(len(unique))

            # Apply the contrast matrix to the identity matrix
            encoded_values = np.dot(identity_matrix, contrast)

            # Create a new DataFrame with the encoded values
            encoded_data = pd.DataFrame(encoded_values, columns=[f"{col}_{i}" for i in unique[:-1]])
            encoded_data[col] = unique
            nonbase_var = col + '_-'
            if nonbase_var in encoded_data.columns:
                encoded_data[nonbase_var] = 0.0
                encoded_data.loc[encoded_data[col] == '-', encoded_data.columns != col] = 0.0
        else:
            encoded_data = pd.DataFrame({col:unique})
        return(encoded_data)

    def get_dropped_contrast_var(self, training_df, original_col):
        """
        Get the names of dropped contrast variables after transformation.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            original_col (str): The name of the original categorical column.

        Returns:
            tuple: A tuple containing two lists - the names of contrast variables and the name of the dropped contrast variable.
        """
        unique = sorted(list(pd.unique(training_df[original_col])))
        contrast_vars=[f"{original_col}_{i}" for i in unique]
        return(contrast_vars[:-1], contrast_vars[-1])

    def transform_categorical_vars(self, training_df, col, pretrain=True):
        """
        Transform categorical variable columns into contrast columns.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            col (str): The name of the categorical column to be transformed.

        Returns:
            pd.DataFrame: The DataFrame with categorical variables transformed into contrast columns.
        """
        assert col in self.variable_colnames, "Input column name is not a variable name"
        first = training_df[col].iloc[0]
        if not isinstance(first, int) and not isinstance(first, float) and not isinstance(first, np.int64):
            contrast_df = self.get_contrast_matrix(training_df, col, pretrain)
            training_df = pd.merge(training_df, contrast_df, on=col, how='inner')
            new_cols = [x for x in list(contrast_df.columns) if x != col]
            self.variable_colnames = [x for x in self.variable_colnames if x != col] + new_cols
            if col in self.choice1_variable_colnames:
                self.choice1_variable_colnames = [x for x in self.choice1_variable_colnames if x != col] + new_cols

            if col in self.choice2_variable_colnames:
                self.choice2_variable_colnames = [x for x in self.choice2_variable_colnames if x != col] + new_cols

        return(training_df)

    def expand_multivariable(self, training_df, col, new_col):
        """
        Expand multi-variable columns into a single string column.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            col (str): The name of the multi-variable column.
            new_col (str): The name of the new string column.

        Returns:
            pd.DataFrame: The DataFrame with the multi-variable column expanded into a string column.
        """
        # Expand multi-variable columns into a single string column
        if type(getattr(self, col)) == list:
            new_col_name = '_'.join(getattr(self, col))
            setattr(self, new_col, new_col_name)
            training_df[new_col_name] = training_df.parallel_apply(lambda row: '_'.join(row[getattr(self, col)].astype(str)), axis=1)
            setattr(self, col, new_col_name)
        return training_df

    def check_within_set_variance(self, training_df, pretrain=True):
        if pretrain:
            for col in self.variable_colnames:
                var_counts = training_df.groupby([self.group_colname, col]).size().reset_index(name = 'N')
                unique_var_counts = var_counts.groupby([self.group_colname]).size().reset_index(name = 'N')
                if 1 in unique_var_counts.N.unique():
                    self.variable_colnames = [item for item in self.variable_colnames if item != col]
                    self.choice1_variable_colnames = [item for item in self.choice1_variable_colnames if item != col]
                    self.choice2_variable_colnames = [item for item in self.choice2_variable_colnames if item != col]
                    print('removing ' + col + ' from model due to insufficient within-choice set variance')
        else:
            # remove missing NT variables
            self.variable_colnames = [var for var in self.variable_colnames if '_-' not in var]
            self.choice1_variable_colnames = [var for var in self.choice1_variable_colnames if '_-' not in var]
            self.choice2_variable_colnames = [var for var in self.choice2_variable_colnames if '_-' not in var]


    def remove_zero_set_counts(self, training_df):
        count_sums = training_df.groupby([self.group_colname])[self.count_colname].sum()
        if 0 in count_sums.unique():
            zeros = count_sums[count_sums == 0]
            zero_groups = list(zeros.index)
            training_df = training_df[~training_df[self.group_colname].isin(zero_groups)]
            print('removing ' + str(zero_groups) + ' groups from model due to zero counts')
        return(training_df)

    def filter_input_domain_space(self, df):
        domain_file = self.params.R_input_domain_data_path()
        domain_data = pd.read_csv(domain_file, sep = '\t')

        # fill in zeros
        # Create a DataFrame of all unique pairs
        unique_genes = domain_data[['v_gene', 'j_gene']].drop_duplicates()
        unique_trims = domain_data[['v_trim', 'j_trim']].drop_duplicates()

        unique_pairs = unique_genes.merge(unique_trims, how='cross')
        unique_pairs['ligation_mh'] = 0

        # merge
        domain_data_subset = domain_data[unique_pairs.columns]
        tog = domain_data_subset.merge(unique_pairs, how='outer').drop_duplicates()

        # filter input df
        filtered_df = pd.merge(tog, df, how='inner', on=tog.columns.tolist())

        return filtered_df


class DataTransformer(DataPreprocessor):
    """
    A class for data preprocessing and transformation.

    Args:
        training_df (pd.DataFrame): The training DataFrame.
        variable_colnames (list): List of variable column names.
        count_colname (str): The column name for the count variable.
        group_colname (str): The column name representing group identifiers.
        repeat_obs_colname (str): The column name representing repeated observation identifiers.
        choice_colname (str): The column name representing choice identifiers.

    Attributes:
        original_variable_colnames (list): List of original variable column names.
        original_group_colname (str): The original column name representing group identifiers.
        original_choice_colname (str): The original column name representing choice identifiers.
        coefs (ndarray): Coefficients for the model.

    Methods:
        get_random_coefs():
            Generates random coefficients for model initialization.

        preprocess_data():
            Preprocesses the training data, transforming columns as necessary.

        get_matrices():
            Prepares data matrices for modeling.

        get_coefficients(coefs=None):
            Retrieves model coefficients as a dictionary.

    Inherits Attributes and Methods from DataPreprocessor class.
    """
    def __init__(self, training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params):
        super().__init__(training_df, variable_colnames, None, None, count_colname, group_colname, repeat_obs_colname, choice_colname, params)
        self.original_variable_colnames = variable_colnames
        self.original_group_colname = group_colname
        self.original_choice_colname = choice_colname
        self.input_variable_colnames = variable_colnames
        self.input_group_colname = group_colname
        self.input_choice_colname = choice_colname
        self.coefs = None

    def get_random_coefs(self):
        """
        Generates random coefficients for model initialization.

        Returns:
            ndarray: Randomly generated coefficients as a NumPy array.
        """
        # Set random coefficients for model initialization
        key = jax.random.PRNGKey(123)
        coefs = jax.random.normal(key, shape=(len(self.variable_colnames), 1))
        return coefs

    def preprocess_data(self, df, pretrain=True):
        """
        Preprocesses the training data, transforming columns as necessary.

        Returns:
            pd.DataFrame: The preprocessed training DataFrame.
        """
        # filter for possible sites
        if 'ligation_mh' in self.input_choice_colname:
            df = self.filter_input_domain_space(df)

        # Transform group column lists into strings
        df = self.expand_multivariable(df, "original_group_colname", "group_colname")

        # remove zero counts
        if pretrain:
            df = self.remove_zero_set_counts(df)

        df = self.transform_categorical_response_vars(df, "original_group_colname", "group_colname")

        # Transform choice column lists into strings
        df = self.expand_multivariable(df, "original_choice_colname", "choice_colname")
        nonrepeat_groups = self.group_colname

        if self.repeat_obs_colname != None:
            self.group_colname = [self.group_colname] + [self.repeat_obs_colname]
            self.original_group_colname = self.group_colname
            df = self.expand_multivariable(df, "original_group_colname", "group_colname")

        # Transform group and choice columns into integer type
        df = self.transform_categorical_response_vars(df, "original_group_colname", "group_colname")
        df = self.transform_categorical_response_vars(df, "original_choice_colname", "choice_colname")

        # Transform categorical variable columns into contrast columns
        for col in self.variable_colnames:
            df = self.transform_categorical_vars(df, col, pretrain)

        # check for within-choice-set variance
        self.check_within_set_variance(df, pretrain)

        return(df)

    def reset_weighted_observations(self, counts_mat):
        mat_sum = jnp.nansum(counts_mat)
        reset_mat = counts_mat/mat_sum
        return(reset_mat)

    def get_matrices(self, df, pretrain=True, replace_object=None, return_df=False):
        """
        Prepares data matrices for modeling.

        Returns:
            tuple: A tuple containing four data matrices - variables, counts, non-repeat groups, and mask.
        """
        df = self.preprocess_data(df, pretrain)

        # Fill in missing counts with zero
        df[self.count_colname] = df[self.count_colname].fillna(0).astype(float)

        # Get matrix shapes
        groups = pd.unique(df[self.group_colname])
        choices = pd.unique(df[self.choice_colname])

        final_shape = (len(groups), len(choices), len(self.variable_colnames))
        int_shape = (len(groups), len(self.variable_colnames), len(choices))
        counts_shape = (len(groups), len(choices), 1)

        # map variables to indices
        var_mapping = {key: value for key, value in enumerate(self.variable_colnames)}
        var_mapping_df = pd.DataFrame(list(var_mapping.items()), columns=['VarIndex', 'VarName'])

        # get counts matrix
        pivot_counts = df.pivot_table(index=[self.group_colname],
                                      columns=[self.choice_colname],
                                      values=self.count_colname,
                                      fill_value=0)
        counts_mat = jnp.array(pivot_counts).reshape(counts_shape)

        # get variable matrix
        pivot_vars = df.pivot_table(index=[self.group_colname],
                                    columns=[self.choice_colname],
                                    values=var_mapping_df.VarName.tolist(),
                                    fill_value=0)
        # reorder columns to reflect correct order
        pivot_vars = pivot_vars.reindex(var_mapping_df.VarName.tolist(), axis=1, level = 0)

        mat = jnp.array(pivot_vars).reshape(int_shape)
        mat = mat.transpose((0, 2, 1))
        assert mat.shape == final_shape, "Variable matrix is the incorrect dimension"

        # Assuming groups is a list of unique group identifiers
        group_indices = {group: idx for idx, group in enumerate(groups)}
        nonrepeat_groups_mat = jnp.array([group_indices[group] for group in df[self.group_colname].unique()])

        # get mask matrix (1 for valid entries or 0 for nonvalid)
        df['indicator'] = 1.0
        pivot_mask = df.pivot_table(index=[self.group_colname],
                                    columns=[self.choice_colname],
                                    values=['indicator'],
                                    fill_value=0)
        mask_mat = jnp.array(pivot_mask).reshape(counts_shape)

        if replace_object is not None:
            setattr(self, replace_object, df)

        if return_df:
            return mat, counts_mat, nonrepeat_groups_mat, mask_mat, df
        else:
            return mat, counts_mat, nonrepeat_groups_mat, mask_mat

    def get_coefficients(self, coefs=None):
        if coefs is None:
            coefs = self.coefs
        coefs = coefs.flatten().tolist()
        errors = self.standard_errors.flatten().tolist()
        choice_sum = self.choice1_variable_colnames + self.choice2_variable_colnames
        if choice_sum != []:
            assert self.variable_colnames == self.choice1_variable_colnames + self.choice2_variable_colnames
        d = {n: {'value':c, 'error':e} for n, c, e in zip(self.variable_colnames, coefs, errors)}
        for col in list(set(self.original_variable_colnames) - set(self.variable_colnames)):
            contrast_vars, missing_var = self.get_dropped_contrast_var(self.training_df, col)
            for element in contrast_vars:
                if element not in self.variable_colnames:
                    contrast_vars.remove(element)
            d[missing_var] = {'value':-1*sum(d[var]['value'] for var in contrast_vars), 'error':None}

        data = [{'coefficient': name, 'value': values['value'], 'error': values['error']} for name, values in d.items()]
        return data

    def get_coefficients_df(self):
        """
        Get the model coefficients as a DataFrame.

        Returns:
            dict: A dictionary containing variable names as keys and their corresponding coefficients as values.
        """
        df = pd.DataFrame(self.get_coefficients())
        df.reset_index(drop=True, inplace=True)

        # get bases
        df['base'] = None
        df.loc[df.coefficient.str.contains('motif'), 'base'] = df.coefficient.str.split('_').str[-1]
        df.loc[df.coefficient.str.contains('base_count'), 'base'] = df.coefficient.str.split('_').str[-1]
        df.loc[df.coefficient.str.contains('base_count') & df.coefficient.str.contains('prop'), 'base'] = df.coefficient.str.split('_').str[-2]


        # get positions
        df['position'] = None
        df.loc[df.coefficient.str.contains('motif'), 'position'] = df.coefficient.str.split('_').str[-2]
        df.loc[df.coefficient.str.contains('mh_prop'), 'position'] = df.coefficient.str.split('_').str[3] + df.coefficient.str.split('_').str[4]
        df.loc[df.coefficient.str.contains('mh_count'), 'position'] = df.coefficient.str.split('_').str[3] + df.coefficient.str.split('_').str[4]


        # get side
        df['side'] = None
        df.loc[df.coefficient.str.contains('motif'), 'side'] = df.coefficient.str.split('_').str[3]
        df.loc[df.coefficient.str.contains('base_count'), 'side'] = df.coefficient.str.split('_').str[2]
        df.loc[df.coefficient.str.contains('mh_prop'), 'side'] = df.coefficient.str.split('_').str[2]
        df.loc[df.coefficient.str.contains('mh_count'), 'side'] = df.coefficient.str.split('_').str[2]

        # fix trimming specific
        df['trim_type'] = None
        df.loc[df.coefficient.str.contains('_trim'), 'trim_type'] = df.coefficient.str.split('_trim').str[0] + '_trim'
        df.loc[df.coefficient.str.contains('_length'), 'trim_type'] = df.coefficient.str.split('_length').str[0].str.split('_').str[-1] + '_trim'

        # simplify coefficients
        df.loc[df.coefficient.str.contains('interaction') & df.coefficient.str.contains('mh_prop'), 'coefficient'] = 'mh_prop_length_interaction'
        df.loc[df.coefficient.str.contains('interaction') & df.coefficient.str.contains('mh_count'), 'coefficient'] = 'mh_count_length_interaction'
        df.loc[df.coefficient.str.contains('motif'), 'coefficient'] = 'motif'
        df.loc[df.coefficient.str.contains('base_count'), 'coefficient'] = 'base_count'
        df.loc[df.coefficient.str.contains('mh_prop') & ~df.coefficient.str.contains('interaction'), 'coefficient'] = 'mh_prop'
        df.loc[df.coefficient.str.contains('mh_count') & ~df.coefficient.str.contains('interaction'), 'coefficient'] = 'mh_count'

        return(df)

class TwoStepDataTransformer(DataTransformer):
    def __init__(self, training_df, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, choice2_colname, params):
        super().__init__(training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params)
        self.choice2_colname = choice2_colname
        self.choice1_variable_colnames = choice1_variable_colnames
        self.choice2_variable_colnames = choice2_variable_colnames
        self.original_choice1_variable_colnames = choice1_variable_colnames
        self.original_choice2_variable_colnames = choice2_variable_colnames
        self.original_choice2_colname = choice2_colname
        self.input_choice1_variable_colnames = choice1_variable_colnames
        self.input_choice2_variable_colnames = choice2_variable_colnames
        self.input_choice2_colname = choice2_colname
        self.coefs = None

    def preprocess_data(self, df, pretrain=True):
        """
        Preprocesses the training data, transforming columns as necessary.

        Returns:
            pd.DataFrame: The preprocessed training DataFrame.
        """
        # filter for possible sites
        if 'ligation_mh' in self.input_choice_colname:
            df = self.filter_input_domain_space(df)
        if 'ligation_mh' in self.input_choice2_colname:
            df = self.filter_input_domain_space(df)

        # Transform group column lists into strings
        df = self.expand_multivariable(df, "original_group_colname", "group_colname")

        # remove zero counts
        if pretrain:
            df = self.remove_zero_set_counts(df)

        df = self.transform_categorical_response_vars(df, "original_group_colname", "group_colname")

        # Transform choice column lists into strings
        df = self.expand_multivariable(df, "original_choice_colname", "choice_colname")
        df = self.expand_multivariable(df, "original_choice2_colname", "choice2_colname")

        nonrepeat_groups = self.group_colname

        if self.repeat_obs_colname != None:
            self.group_colname = [self.group_colname] + [self.repeat_obs_colname]
            self.original_group_colname = self.group_colname
            df = self.expand_multivariable(df, "original_group_colname", "group_colname")

        # Transform group and choice columns into integer type
        df = self.transform_categorical_response_vars(df, "original_group_colname", "group_colname")
        df = self.transform_categorical_response_vars(df, "original_choice_colname", "choice_colname")
        df = self.transform_categorical_response_vars(df, "original_choice2_colname", "choice2_colname")


        # Transform categorical variable columns into contrast columns
        for col in self.variable_colnames:
            df = self.transform_categorical_vars(df, col, pretrain)

        # check for within-choice-set variance
        self.check_within_set_variance(df, pretrain)

        return(df)


    def get_matrices(self, df, pretrain=True, replace_object=None, return_df=False):
        """
        Prepares data matrices for twostep modeling.

        Returns:
            tuple: A tuple containing four data matrices - variables, counts, non-repeat groups, and mask.
        """
        #TODO may need to create a second mask for nonproductive filtering!!!
        df = self.preprocess_data(df, pretrain)

        # Fill in missing counts with zero
        df[self.count_colname] = df[self.count_colname].fillna(0).astype(float)

        # Get matrix shapes
        groups = pd.unique(df[self.group_colname])
        choices = pd.unique(df[self.choice_colname])
        choice2s = pd.unique(df[self.choice2_colname])

        choice1_final_shape = (len(groups), len(choices), len(self.choice1_variable_colnames))
        choice1_int_shape = (len(groups), len(self.choice1_variable_colnames), len(choices))
        choice2_final_shape = (len(groups), len(choices), len(choice2s), len(self.choice2_variable_colnames))
        choice2_int_shape = (len(groups), len(self.choice2_variable_colnames), len(choices), len(choice2s))
        counts_shape = (len(groups), len(choices), len(choice2s), 1)

        # map variables to indices
        choice1_var_mapping = {key: value for key, value in enumerate(self.choice1_variable_colnames)}
        choice1_var_mapping_df = pd.DataFrame(list(choice1_var_mapping.items()), columns=['VarIndex', 'VarName'])
        choice2_var_mapping = {key: value for key, value in enumerate(self.choice2_variable_colnames)}
        choice2_var_mapping_df = pd.DataFrame(list(choice2_var_mapping.items()), columns=['VarIndex', 'VarName'])


        # fill in df with all combos
        df['indicator'] = 1.0
        all_combinations = pd.MultiIndex.from_product([choices, choice2s],
                                                      names=[self.choice_colname, self.choice2_colname]).to_frame(index=False)
        df_full = pd.merge(all_combinations, df, on=[self.choice_colname, self.choice2_colname], how='left').fillna(0)

        # get counts matrix
        pivot_counts = df_full.pivot_table(index=[self.group_colname],
                                      columns=[self.choice_colname, self.choice2_colname],
                                      values=self.count_colname,
                                      fill_value=0)
        counts_mat = jnp.array(pivot_counts).reshape(counts_shape)

        # get variable matrix
        choice1_pivot_vars = df_full.pivot_table(index=[self.group_colname],
                                    columns=[self.choice_colname],
                                    values=choice1_var_mapping_df.VarName.tolist(),
                                    fill_value=0)
        # reorder columns to reflect correct order
        choice1_pivot_vars = choice1_pivot_vars.reindex(choice1_var_mapping_df.VarName.tolist(), axis=1, level = 0)

        choice2_pivot_vars = df_full.pivot_table(index=[self.group_colname],
                                    columns=[self.choice_colname, self.choice2_colname],
                                    values=choice2_var_mapping_df.VarName.tolist(),
                                    fill_value=0)
        # reorder columns to reflect correct order
        choice2_pivot_vars = choice2_pivot_vars.reindex(choice2_var_mapping_df.VarName.tolist(), axis=1, level = 0)

        choice1_mat = jnp.array(choice1_pivot_vars).reshape(choice1_int_shape)
        choice1_mat = choice1_mat.transpose((0, 2, 1))
        assert choice1_mat.shape == choice1_final_shape, "choice1 variable matrix is the incorrect dimension"

        choice2_mat = jnp.array(choice2_pivot_vars).reshape(choice2_int_shape)
        choice2_mat = choice2_mat.transpose((0, 2, 3, 1))
        assert choice2_mat.shape == choice2_final_shape, "choice2 variable matrix is the incorrect dimension"

        # Assuming groups is a list of unique group identifiers
        group_indices = {group: idx for idx, group in enumerate(groups)}
        nonrepeat_groups_mat = jnp.array([group_indices[group] for group in df[self.group_colname].unique()])

        # get mask matrix (1 for valid entries or 0 for nonvalid)
        pivot_mask = df_full.pivot_table(index=[self.group_colname],
                                    columns=[self.choice_colname, self.choice2_colname],
                                    values=['indicator'],
                                    fill_value=0)
        mask_mat = jnp.array(pivot_mask).reshape(counts_shape)

        if replace_object is not None:
            setattr(self, replace_object, df)

        if return_df:
            return choice1_mat, choice2_mat, counts_mat, nonrepeat_groups_mat, mask_mat, df
        else:
            return choice1_mat, choice2_mat, counts_mat, nonrepeat_groups_mat, mask_mat

    def get_random_coefs(self):
        """
        Generates random coefficients for model initialization.

        Returns:
            ndarray: Randomly generated coefficients as a NumPy array.
        """
        # Set random coefficients for model initialization
        key = jax.random.PRNGKey(123)
        coefs = jax.random.normal(key, shape=(len(self.choice1_variable_colnames) + len(self.choice2_variable_colnames), 1))
        return coefs


class ConditionalLogisticRegressor(DataTransformer):
    """
    A class for training a conditional logistic regression model and performing cross-validation for hyperparameter tuning.

    Args:
        training_df (pd.DataFrame): The training DataFrame.
        variable_colnames (list): List of variable column names.
        count_colname (str): The column name for the count variable.
        group_colname (str): The column name representing group identifiers.
        repeat_obs_colname (str): The column name representing repeated observation identifiers.
        choice_colname (str): The column name representing choice identifiers.
        l2kfold (int): Number of folds for L2 regularization hyperparameter tuning.

    Attributes:
        original_variable_colnames (list): List of original variable column names.
        original_group_colname (str): The original column name representing group identifiers.
        original_choice_colname (str): The original column name representing choice identifiers.
        variable_matrix (ndarray): Data matrix of variables.
        counts_matrix (ndarray): Data matrix of counts.
        nonrepeat_grp_matrix (ndarray): Data matrix of non-repeating groups.
        initial_coefs (ndarray): Initial coefficients for model training.
        coefs (ndarray): Trained model coefficients.
        training_info: Information about the training process.
        maxiter (int): Maximum number of iterations for optimization.
        tolerance (float): tolerance for optimization.
        l2reg (float): L2 regularization strength.
        l2kfold (int): Number of folds for L2 regularization hyperparameter tuning.

    Methods:
        get_prob(variables, mask, coefs=None):
            Computes choice probabilities for given variables and coefficients.

        cross_entropy(probs, counts):
            Calculates the cross-entropy loss.

        l2regularization(coefs, size, l2reg):
            Computes L2 regularization term.

        loss_fn(coefs, variables, counts, mask, l2reg=0):
            Computes the loss function for optimization.

        fit(variable_matrix, counts_matrix, mask_matrix, l2reg, maxiter, tolerance, initial_coefs):
            Fits the conditional logistic regression model using gradient descent.

        grid_search_cv(l2kfold, l2reg_values=[10**i for i in range(-2, 8)] + [0]):
            Performs grid search cross-validation for hyperparameter tuning.

        get_l2reg():
            Determines the optimal L2 regularization strength.

        train_model(l2=False, maxiter=None, tolerance=None):
            Trains the conditional logistic regression model with optional L2 regularization.

    Inherits Attributes and Methods from DataTransformer class.
    """
    def __init__(self, training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params, l2kfold=10):
        super().__init__(training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params)
        self.variable_matrix, self.counts_matrix, self.nonrepeat_grp_matrix, self.mask_matrix = self.get_matrices(training_df, replace_object='training_df')
        self.initial_coefs = self.get_random_coefs()
        self.coefs = None
        self.training_info = None
        self.maxiter = 1000
        self.tolerance = 1e-6
        self.step = 0.1
        self.l2reg = 0
        self.l2kfold = None
        self.l2reg_grid = None

    # Get probability for input parameters given coefficients
    def get_prob(self, variables, mask, coefs=None):
        """
        Compute choice probabilities for given variables and coefficients.

        Args:
            variables (ndarray): Data matrix of variables.
            mask (ndarray): Data matrix of variable mask.
            coefs (ndarray, optional): Coefficients for the logistic regression model. If not provided, the trained coefficients will be used.

        Returns:
            ndarray: Choice probabilities for each choice given the variables.
        """
        if coefs is None:
            coefs = self.coefs
        assert coefs is not None, "Need to train the model before making predictions!"
        # Compute the logits for each choice
        cov = jnp.dot(variables, coefs)
        reshape = jnp.squeeze(cov)
        reshaped_mask = jnp.squeeze(mask)

        # Calculate the probability of the observed choices
        # Dimensions of this matrix are groups x choices
        # replace missing choices with -INF so that they will not count towards probability
        probs = jax.nn.softmax(jnp.where(reshaped_mask, reshape, jnp.NINF))

        return probs

    # Get cross-entropy loss
    def cross_entropy(self, probs, counts):
        """
        Calculate the cross-entropy loss.

        Args:
            probs (ndarray): Choice probabilities.
            counts (ndarray): Counts of choices.

        Returns:
            float: Cross-entropy loss.
        """
        counts_reshape = jnp.squeeze(counts)
        loss = -jnp.sum(jnp.log(jnp.where(probs==0, 1, probs)) * counts_reshape)
        return loss

    def l2regularization(self, coefs, size, l2reg):
        """
        Compute L2 regularization term for only non-base-count and non-motif coefficients.

        Args:
            coefs (ndarray): Coefficients for the logistic regression model.
            size (int): Size of the training data
            l2reg (float): L2 regularization strength.

        Returns:
            float: L2 regularization term.
        """
        var_list = [var for var in self.variable_colnames if 'base_count' not in var or 'interaction' in var]
        var_list = [var for var in var_list if 'motif' not in var]

        mh_list = [var in var_list for var in self.variable_colnames]

        if any(mh_list):
            binary_mh_list = [int(b) for b in mh_list]
            binary_mh_jnp = jnp.array(binary_mh_list).reshape(-1, 1)
            coef_subset = coefs * binary_mh_jnp
            c = jnp.nansum(coef_subset**2)
        else:
            c = 0
        return(0.5*(1/size)*l2reg*c)

    # Compute the loss function
    def loss_fn(self, coefs, variables, counts, mask, l2reg=0):
        """
        Compute the loss function for optimization.

        Args:
            coefs (ndarray): Coefficients for the logistic regression model.
            variables (ndarray): Data matrix of variables.
            counts (ndarray): Counts of choices.
            mask (ndarray): Data matrix of mask.
            l2reg (float): L2 regularization strength.

        Returns:
            float: Total loss, including cross-entropy and L2 regularization.
        """
        probs = self.get_prob(variables, mask, coefs)
        if probs.ndim == 1:
            probs = probs.reshape((probs.shape[0], 1))
        size = jnp.count_nonzero(mask).item()

        loss = self.cross_entropy(probs, counts) + self.l2regularization(coefs, size, l2reg)
        return loss

    def fit(self, variable_matrix, counts_matrix, mask_matrix, l2reg, maxiter, tol, step, initial_coefs):
        """
        Fit the conditional logistic regression model using gradient descent.

        Args:
            variable_matrix (ndarray): Data matrix of variables.
            counts_matrix (ndarray): Counts of choices.
            mask_matrix (ndarray): Data matrix of masks.
            l2reg (float): L2 regularization strength.
            maxiter (int): Maximum number of optimization iterations.
            tolerance (float): tolerance for optimization.
            initial_coefs (ndarray): Initial coefficients for model training.

        Returns:
            OptimizationResult: Result of the optimization process.
        """
        assert counts_matrix is not None, "counts column is missing"

        # Create a jaxopt GradientDescent optimizer
        solver = jaxopt.BFGS(fun=self.loss_fn, maxiter=maxiter, tol=tol, verbose=True)

        # Run gradient descent
        res = solver.run(initial_coefs,
                         variables=variable_matrix,
                         counts=counts_matrix,
                         mask=mask_matrix,
                         l2reg=l2reg)
        return(res)

    def cv_loss(self, fold_count, l2reg):
        assert self.counts_matrix is not None, "counts column is needed"

        kf = GroupKFold(n_splits=fold_count)
        scores = []

        for train_index, val_index in kf.split(X=self.variable_matrix, y=self.counts_matrix, groups=self.nonrepeat_grp_matrix):
            train_data, val_data = self.variable_matrix[train_index], self.variable_matrix[val_index]
            train_counts, val_counts = self.counts_matrix[train_index], self.counts_matrix[val_index]
            train_mask, val_mask = self.mask_matrix[train_index], self.mask_matrix[val_index]
            train_counts = self.reset_weighted_observations(train_counts)
            val_counts = self.reset_weighted_observations(val_counts)

            # Train the model on the training data
            model = self.fit(train_data, train_counts, train_mask, l2reg, self.maxiter, self.tolerance, self.step, self.initial_coefs)

            # Compute the loss on the validation data
            loss = self.loss_fn(model.params,
                                val_data,
                                val_counts,
                                val_mask)

            # Store the loss as a score (lower score is better)
            scores.append(float(loss))
        return(scores)

    def grid_search_cv(self, l2kfold, l2reg_values=[10**i for i in range(-5, 5)] + [0]):
        """
        Perform grid search cross-validation for hyperparameter tuning.

        Args:
            l2kfold (int): Number of folds for cross-validation.
            l2reg_values (list): List of L2 regularization strengths to search.

        Returns:
            pd.DataFrame: DataFrame with hyperparameter values and corresponding mean cross-validation loss.
        """
        results = pd.DataFrame({'l2reg':[], 'mean_CV_loss':[]})

        for l2reg in l2reg_values:
            # Perform cross-validation
            scores = self.cv_loss(l2kfold, l2reg)

            # Calculate the mean cross-validation score
            mean_score = np.mean(scores)

            temp = pd.DataFrame({'l2reg':[l2reg], 'mean_CV_loss':[mean_score]})
            results = pd.concat([results, temp], ignore_index=True)

        return(results)

    def get_l2reg(self):
        """
        Determine the optimal L2 regularization strength.

        Returns:
            float: Optimal L2 regularization strength.
        """
        self.l2kfold = min(len(np.unique(self.nonrepeat_grp_matrix)), 10)
        grid = self.grid_search_cv(self.l2kfold)
        self.l2reg_grid = grid
        min_obs = grid[grid.mean_CV_loss == min(grid.mean_CV_loss)]
        return(float(min_obs.l2reg.iloc[0]))

    def train_model(self, l2=False, l2reg_value=None, maxiter=None, tolerance=None, step=None):
        """
        Train the conditional logistic regression model with optional L2 regularization.

        Args:
            l2 (bool): Whether to use L2 regularization.
            maxiter (int, optional): Maximum number of optimization iterations. If not provided, the default value is used.
            tolerance (float, optional): tolerance for optimization. If not provided, the default value is used.

        Returns:
            self: The trained ConditionalLogisticRegressor instance.
        """
        assert self.counts_matrix is not None, "counts column is needed"

        if l2 is False:
            self.l2reg = 0
        else:
            if l2reg_value is None:
                self.l2reg = self.get_l2reg()
            else:
                self.l2reg = l2reg_value

        if maxiter is None:
            maxiter = self.maxiter
        else:
            self.maxiter = maxiter

        if tolerance is None:
            tolerance = self.tolerance
        else:
            self.tolerance = tolerance

        if step is None:
            step = self.step
        else:
            self.step = step

        res = self.fit(self.variable_matrix,
                       self.counts_matrix,
                       self.mask_matrix,
                       self.l2reg,
                       maxiter,
                       tolerance,
                       step,
                       self.initial_coefs)

        self.coefs = res.params
        self.training_info = res
        self.maxiter = maxiter
        self.tolerance = tolerance
        self.step = step
        self.cov_matrix = self.get_cov_matrix(self.coefs, self.variable_matrix, self.counts_matrix, self.mask_matrix, self.l2reg)
        self.standard_errors = self.get_errors(self.coefs, self.variable_matrix, self.counts_matrix, self.mask_matrix, self.l2reg)
        return self

    def get_hessian(self, coefs, variables, counts, mask, l2reg=0):
         # Wrapper function
        def wrapper_loss_fn(coefs):
            return self.loss_fn(coefs, variables, counts, mask, l2reg)

        hessian_fn = jax.hessian(wrapper_loss_fn, argnums=0)
        hessian_matrix = hessian_fn(coefs.reshape(-1))
        return hessian_matrix

    def get_cov_matrix(self, coefs, variables, counts, mask, l2reg=0):
        hess_mat = self.get_hessian(coefs, variables, counts, mask, l2reg)
        cov_matrix = jnp.linalg.inv(hess_mat)
        return cov_matrix

    def get_errors(self, coefs, variables, counts, mask, l2reg=0):
        cov = self.get_cov_matrix(coefs, variables, counts, mask, l2reg)
        standard_errors = np.sqrt(np.diag(cov))
        return standard_errors

    def save_model(self, file_path):
        assert self.coefs is not None, "need to train model before saving"
        self.training_df = None
        self.variable_matrix = None
        self.counts_matrix = None
        self.nonrepeat_grp_matrix = None
        with open(file_path, 'wb') as file:
            # Serialize and save the object to the file
            dill.dump(self, file)

class TwoStepConditionalLogisticRegressor(TwoStepDataTransformer):
    def __init__(self, training_df, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, choice2_colname, params, l2kfold=10):
        super().__init__(training_df, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, choice2_colname, params)
        self.choice1_variable_matrix, self.choice2_variable_matrix, self.counts_matrix, self.nonrepeat_grp_matrix, self.mask_matrix = self.get_matrices(training_df, replace_object='training_df')
        self.initial_coefs = self.get_random_coefs()
        self.coefs = None
        self.training_info = None
        self.maxiter = 1000
        self.tolerance = 1e-6
        self.step = 0.1
        self.l2reg = 0
        self.l2kfold = None
        self.l2reg_grid = None

    # Get probability for input parameters given coefficients
    def get_indiv_prob(self, choice1_variables, choice2_variables, mask, coefs=None):
        if coefs is None:
            coefs = self.coefs

        assert coefs is not None, "Need to train the model before making predictions!"

        choice1_coefs = coefs[0:len(self.choice1_variable_colnames)]
        choice2_coefs = coefs[len(self.choice1_variable_colnames):]

        # Compute the logits for each choice
        choice1_cov = jnp.dot(choice1_variables, choice1_coefs)
        choice1_reshape = jnp.squeeze(choice1_cov)
        choice2_cov = jnp.dot(choice2_variables, choice2_coefs)
        choice2_reshape = jnp.squeeze(choice2_cov)
        reshaped_mask = jnp.squeeze(mask)

        # Calculate the probability of the observed choices
        # Dimensions of this matrix are groups x choices
        # replace missing choices with -INF so that they will not count towards probability
        choice1_probs = jax.nn.softmax(choice1_reshape)
        choice2_probs = jax.nn.softmax(jnp.where(reshaped_mask, choice2_reshape, jnp.NINF))

        return choice1_probs, choice2_probs

    def get_joint_prob(self, choice1_variables, choice2_variables, mask, coefs=None):
        choice1_probs, choice2_probs = self.get_indiv_prob(choice1_variables, choice2_variables, mask, coefs)
        choice1_probs_reshape = choice1_probs.reshape(choice1_probs.shape[0], choice1_probs.shape[1], 1)
        prob = choice1_probs_reshape * choice2_probs
        return prob

    # Get cross-entropy loss
    def cross_entropy(self, choice1_probs, choice2_probs, counts):
        """
        Calculate the cross-entropy loss.

        Args:
            probs (ndarray): Choice probabilities.
            counts (ndarray): Counts of choices.

        Returns:
            float: Cross-entropy loss.
        """
        counts1 = counts.sum(axis = 2)
        counts2 = counts
        counts1_reshape = jnp.squeeze(counts1)
        counts2_reshape = jnp.squeeze(counts2)

        loss1 = -jnp.sum(jnp.log(jnp.where(choice1_probs==0, 1, choice1_probs)) * counts1_reshape)
        loss2 = -jnp.sum(jnp.log(jnp.where(choice2_probs==0, 1, choice2_probs)) * counts2_reshape)

        loss = loss1 +  loss2

        return loss

    def l2regularization(self, coefs, size, l2reg):
        """
        Compute L2 regularization term for only non-base-count and non-motif coefficients.

        Args:
            coefs (ndarray): Coefficients for the logistic regression model.
            size (int): Size of the training data
            l2reg (float): L2 regularization strength.

        Returns:
            float: L2 regularization term.
        """
        def calculate_coef_sum(choice_var_colnames, coefs):
            var_list = [var for var in choice_var_colnames if 'base_count' not in var or 'interaction' in var]
            var_list = [var for var in var_list if 'motif' not in var]

            mh_list = [var in var_list for var in choice_var_colnames]
            mh_jnp = jnp.array(mh_list).reshape(-1, 1)

            coef_subset = coefs * mh_jnp

            c = jnp.nansum(coef_subset**2)
            return c

        c = 0
        # add coefs from choice1 variables
        c += calculate_coef_sum(self.choice1_variable_colnames + self.choice2_variable_colnames, coefs)
        return(0.5*(1/size)*l2reg*c)

    # Compute the loss function
    def loss_fn(self, coefs, choice1_variables, choice2_variables, counts, mask, l2reg=0):
        """
        Compute the loss function for optimization.

        Args:
            coefs (ndarray): Coefficients for the logistic regression model.
            variables (ndarray): Data matrix of variables.
            counts (ndarray): Counts of choices.
            mask (ndarray): Data matrix of mask.
            l2reg (float): L2 regularization strength.

        Returns:
            float: Total loss, including cross-entropy and L2 regularization.
        """
        choice1_probs, choice2_probs = self.get_indiv_prob(choice1_variables, choice2_variables, mask, coefs)
        if choice1_probs.ndim == 1:
            choice1_probs = choice1_probs.reshape((choice1_probs.shape[0], 1))
        if choice2_probs.ndim == 1:
            choice2_probs = choice2_probs.reshape((choice2_probs.shape[0], 1))

        size = jnp.count_nonzero(mask).item()

        loss = self.cross_entropy(choice1_probs, choice2_probs, counts) + self.l2regularization(coefs, size, l2reg)
        return loss

    def fit(self, choice1_variable_matrix, choice2_variable_matrix, counts_matrix, mask_matrix, l2reg, maxiter, tol, step, initial_coefs):
        """
        Fit the conditional logistic regression model using gradient descent.

        Args:
            variable_matrix (ndarray): Data matrix of variables.
            counts_matrix (ndarray): Counts of choices.
            mask_matrix (ndarray): Data matrix of masks.
            l2reg (float): L2 regularization strength.
            maxiter (int): Maximum number of optimization iterations.
            tolerance (float): tolerance for optimization.
            initial_coefs (ndarray): Initial coefficients for model training.

        Returns:
            OptimizationResult: Result of the optimization process.
        """
        assert counts_matrix is not None, "counts column is missing"

        # Create a jaxopt GradientDescent optimizer
        solver = jaxopt.BFGS(fun=self.loss_fn, maxiter=maxiter, tol=tol, verbose=True)

        # Run gradient descent
        res = solver.run(initial_coefs,
                         choice1_variables=choice1_variable_matrix,
                         choice2_variables=choice2_variable_matrix,
                         counts=counts_matrix,
                         mask=mask_matrix,
                         l2reg=l2reg)
        return(res)

    def cv_loss(self, fold_count, l2reg):
        assert self.counts_matrix is not None, "counts column is needed"

        kf = GroupKFold(n_splits=fold_count)
        scores = []

        for train_index, val_index in kf.split(X=self.choice1_variable_matrix, y=self.counts_matrix, groups=self.nonrepeat_grp_matrix):
            train_c1_data, train_c2_data, val_c1_data, val_c2_data = self.choice1_variable_matrix[train_index], self.choice2_variable_matrix[train_index], self.choice1_variable_matrix[val_index], self.choice2_variable_matrix[val_index]
            train_counts, val_counts = self.counts_matrix[train_index], self.counts_matrix[val_index]
            train_mask, val_mask = self.mask_matrix[train_index], self.mask_matrix[val_index]
            train_counts = self.reset_weighted_observations(train_counts)
            val_counts = self.reset_weighted_observations(val_counts)

            # Train the model on the training data
            model = self.fit(train_c1_data, train_c2_data, train_counts, train_mask, l2reg, self.maxiter, self.tolerance, self.step, self.initial_coefs)

            # Compute the loss on the validation data
            loss = self.loss_fn(model.params,
                                val_c1_data,
                                val_c2_data,
                                val_counts,
                                val_mask)

            # Store the loss as a score (lower score is better)
            scores.append(float(loss))
        return(scores)

    def grid_search_cv(self, l2kfold, l2reg_values=[10**i for i in range(-5, 5)] + [0]):
        """
        Perform grid search cross-validation for hyperparameter tuning.

        Args:
            l2kfold (int): Number of folds for cross-validation.
            l2reg_values (list): List of L2 regularization strengths to search.

        Returns:
            pd.DataFrame: DataFrame with hyperparameter values and corresponding mean cross-validation loss.
        """
        results = pd.DataFrame({'l2reg':[], 'mean_CV_loss':[]})

        for l2reg in l2reg_values:
            # Perform cross-validation
            scores = self.cv_loss(l2kfold, l2reg)

            # Calculate the mean cross-validation score
            mean_score = np.mean(scores)

            temp = pd.DataFrame({'l2reg':[l2reg], 'mean_CV_loss':[mean_score]})
            results = pd.concat([results, temp], ignore_index=True)

        return(results)

    def get_l2reg(self):
        """
        Determine the optimal L2 regularization strength.

        Returns:
            float: Optimal L2 regularization strength.
        """
        self.l2kfold = min(len(np.unique(self.nonrepeat_grp_matrix)), 10)
        grid = self.grid_search_cv(self.l2kfold)
        self.l2reg_grid = grid
        min_obs = grid[grid.mean_CV_loss == min(grid.mean_CV_loss)]
        return(float(min_obs.l2reg.iloc[0]))

    def train_model(self, l2=False, l2reg_value=None, maxiter=None, tolerance=None, step=None):
        """
        Train the conditional logistic regression model with optional L2 regularization.

        Args:
            l2 (bool): Whether to use L2 regularization.
            maxiter (int, optional): Maximum number of optimization iterations. If not provided, the default value is used.
            tolerance (float, optional): tolerance for optimization. If not provided, the default value is used.

        Returns:
            self: The trained ConditionalLogisticRegressor instance.
        """
        assert self.counts_matrix is not None, "counts column is needed"

        if l2 is False:
            self.l2reg = 0
        else:
            if l2reg_value is None:
                self.l2reg = self.get_l2reg()
            else:
                self.l2reg = l2reg_value

        if maxiter is None:
            maxiter = self.maxiter
        else:
            self.maxiter = maxiter

        if tolerance is None:
            tolerance = self.tolerance
        else:
            self.tolerance = tolerance

        if step is None:
            step = self.step
        else:
            self.step = step

        res = self.fit(self.choice1_variable_matrix,
                       self.choice2_variable_matrix,
                       self.counts_matrix,
                       self.mask_matrix,
                       self.l2reg,
                       maxiter,
                       tolerance,
                       step,
                       self.initial_coefs)

        self.coefs = res.params
        self.training_info = res
        self.maxiter = maxiter
        self.tolerance = tolerance
        self.step = step
        self.cov_matrix = self.get_cov_matrix(self.coefs, self.choice1_variable_matrix, self.choice2_variable_matrix, self.counts_matrix, self.mask_matrix, self.l2reg)
        self.standard_errors = self.get_errors(self.coefs, self.choice1_variable_matrix, self.choice2_variable_matrix, self.counts_matrix, self.mask_matrix, self.l2reg)
        return self

    def get_hessian(self, coefs, choice1_variables, choice2_variables, counts, mask, l2reg=0):
         # Wrapper function
        def wrapper_loss_fn(coefs):
            return self.loss_fn(coefs, choice1_variables, choice2_variables, counts, mask, l2reg)

        hessian_fn = jax.hessian(wrapper_loss_fn, argnums=0)
        hessian_matrix = hessian_fn(coefs.reshape(-1))
        return hessian_matrix

    def get_cov_matrix(self, coefs, choice1_variables, choice2_variables, counts, mask, l2reg=0):
        hess_mat = self.get_hessian(coefs, choice1_variables, choice2_variables, counts, mask, l2reg)
        cov_matrix = jnp.linalg.inv(hess_mat)
        return cov_matrix

    def get_errors(self, coefs, choice1_variables, choice2_variables, counts, mask, l2reg=0):
        cov = self.get_cov_matrix(coefs, choice1_variables, choice2_variables, counts, mask, l2reg)
        standard_errors = np.sqrt(np.diag(cov))
        return standard_errors

    def save_model(self, file_path):
        assert self.coefs is not None, "need to train model before saving"
        self.training_df = None
        self.choice1_variable_matrix = None
        self.choice2_variable_matrix = None
        self.counts_matrix = None
        self.nonrepeat_grp_matrix = None
        with open(file_path, 'wb') as file:
            # Serialize and save the object to the file
            dill.dump(self, file)


class ConditionalLogisticRegressionPredictor(DataTransformer):
    """
    A class for making predictions using a trained conditional logistic regression model.

    Args:
        model (ConditionalLogisticRegressor): A trained conditional logistic regression model.
        variable_colnames (list): List of variable column names.
        count_colname (str): The column name for the count variable.
        group_colname (list): List of column names representing group identifiers.
        choice_colname (list): List of column names representing choice identifiers.

    Attributes:
        model (ConditionalLogisticRegressor): The trained conditional logistic regression model.
        original_variable_colnames (list): List of original variable column names.
        original_group_colname (list): List of original column names representing group identifiers.
        original_choice_colname (list): List of original column names representing choice identifiers.
        variable_matrix (numpy.ndarray): The matrix of predictor variables for prediction.
        counts_matrix (numpy.ndarray): The matrix of counts for prediction.
        nonrepeat_grp_matrix (numpy.ndarray): The matrix of non-repeated group identifiers.

    Methods:
        predict(): Make predictions using the model.
        compute_loss(): Compute the loss for the given data.
        get_coefficients(): Get the model coefficients as a dictionary.
    """
    def __init__(self, model, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params):
        super().__init__(None, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params)
        self.model = model
        if not isinstance(model, ConditionalLogisticRegressor):
            raise TypeError("'model' must be a ConditionalLogisticRegressor object")
        self.original_variable_colnames = variable_colnames
        self.original_group_colname = group_colname
        self.original_choice_colname = choice_colname
        self.input_variable_colnames = variable_colnames
        self.input_group_colname = group_colname
        self.input_choice_colname = choice_colname

    # get probability for input parameters given coefficients
    def predict(self, new_df):
        """
        Make predictions using the trained model.

        Returns:
            pd.DataFrame: A DataFrame containing predicted probabilities merged with the original data.
        Raises:
            ValueError: If the number of variable columns in the input DataFrame doesn't match the model's coefficients.
        """
        original_new_df = new_df
        variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix, new_df = self.get_matrices(new_df, pretrain=False, return_df=True)

        if not variable_matrix.shape[-1] == self.model.coefs.shape[0]:
            raise ValueError("Input dataframe variable column count doesn't match the trained model coefficient count")
        # get predicted probabilities
        probs = self.model.get_prob(variable_matrix, mask_matrix, self.model.coefs)
        # transform probs to a dataframe
        choice_cols = self.get_mapping_dict(new_df, self.choice_colname)
        group_cols = self.get_mapping_dict(new_df, self.group_colname)
        if len(group_cols) == 1:
            probs = probs.reshape((len(group_cols), len(choice_cols)))
        prob_df = pd.DataFrame(probs)
        prob_df.columns = list(choice_cols.keys())
        prob_df[self.group_colname] = list(group_cols.keys())
        melted_df = pd.melt(prob_df,
                            id_vars=[self.group_colname],
                            var_name=self.choice_colname,
                            value_name='predicted_prob')
        # merge predicted probabilities with original df
        merged_df = pd.merge(new_df, melted_df,
                             on=[self.group_colname, self.choice_colname],
                             how='inner')

        # fill in out of domain values
        common_cols = list(set(original_new_df) & set(merged_df))
        final_df = pd.merge(merged_df, original_new_df, how='outer', on=common_cols)

        # Fill missing values with a specified value
        final_df.fillna(0, inplace=True)
        return(final_df)

    def compute_loss(self, new_df):
        """
        Compute the loss for the given data.

        Returns:
            float: The computed loss.
        """
        variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(new_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        loss = self.model.loss_fn(self.model.coefs,
                                  variable_matrix,
                                  counts_matrix,
                                  mask_matrix)
        return(float(loss))

class TwoStepConditionalLogisticRegressionPredictor(TwoStepDataTransformer):
    """
    A class for making predictions using a trained conditional logistic regression model.

    Args:
        model (ConditionalLogisticRegressor): A trained conditional logistic regression model.
        variable_colnames (list): List of variable column names.
        count_colname (str): The column name for the count variable.
        group_colname (list): List of column names representing group identifiers.
        choice_colname (list): List of column names representing choice identifiers.

    Attributes:
        model (ConditionalLogisticRegressor): The trained conditional logistic regression model.
        original_variable_colnames (list): List of original variable column names.
        original_group_colname (list): List of original column names representing group identifiers.
        original_choice_colname (list): List of original column names representing choice identifiers.
        variable_matrix (numpy.ndarray): The matrix of predictor variables for prediction.
        counts_matrix (numpy.ndarray): The matrix of counts for prediction.
        nonrepeat_grp_matrix (numpy.ndarray): The matrix of non-repeated group identifiers.

    Methods:
        predict(): Make predictions using the model.
        compute_loss(): Compute the loss for the given data.
        get_coefficients(): Get the model coefficients as a dictionary.
    """
    def __init__(self, model, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, choice2_colname, params):
        super().__init__(None, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, choice2_colname, params)
        self.model = model
        if not isinstance(model, TwoStepConditionalLogisticRegressor):
            raise TypeError("'model' must be a TwoStepConditionalLogisticRegressor object")

    # get probability for input parameters given coefficients
    def predict(self, new_df):
        """
        Make predictions using the trained model.

        Returns:
            pd.DataFrame: A DataFrame containing predicted probabilities merged with the original data.
        Raises:
            ValueError: If the number of variable columns in the input DataFrame doesn't match the model's coefficients.
        """
        original_new_df = new_df
        choice1_variable_matrix, choice2_variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix, new_df = self.get_matrices(new_df, pretrain=False, return_df=True)

        if not choice1_variable_matrix.shape[-1] + choice2_variable_matrix.shape[-1] == self.model.coefs.shape[0]:
            raise ValueError("Input dataframe variable column count doesn't match the trained model coefficient count")
        # get predicted probabilities
        probs = self.model.get_joint_prob(choice1_variable_matrix, choice2_variable_matrix, mask_matrix, self.model.coefs)
        # transform probs to a dataframe
        choice1_cols = self.get_mapping_dict(new_df, self.choice_colname)
        choice2_cols = self.get_mapping_dict(new_df, self.choice2_colname)
        group_cols = self.get_mapping_dict(new_df, self.group_colname)
        if len(group_cols) == 1:
            probs = probs.reshape((len(group_cols), len(choice1_cols), len(choice2_cols)))

        # Convert the 3D array to a 2D array where each "page" is flattened out
        # and then stack them vertically
        array_2d = probs.reshape(-1, probs.shape[2])

        # Create a MultiIndex representing the first two dimensions (depth and rows)
        # np.repeat and np.tile help to properly align the indices with the reshaped array
        index = pd.MultiIndex.from_tuples([(i, j) for i in range(probs.shape[0]) for j in range(probs.shape[1])], names=[self.group_colname, self.choice_colname])

        # Create the DataFrame
        prob_df = pd.DataFrame(array_2d, index=index, columns=[i for i in range(probs.shape[2])])
        prob_df = prob_df.reset_index()

        melted_df = pd.melt(prob_df,
                            id_vars=[self.group_colname, self.choice_colname],
                            var_name=self.choice2_colname,
                            value_name='predicted_prob')

        # Invert the dictionary
        choice1_mapping = {v: k for k, v in choice1_cols.items()}
        choice2_mapping = {v: k for k, v in choice2_cols.items()}
        group_mapping = {v: k for k, v in group_cols.items()}

        # Use the map function to replace values in 'Column1'
        melted_df[self.group_colname] = melted_df[self.group_colname].map(group_mapping)
        melted_df[self.choice_colname] = melted_df[self.choice_colname].map(choice1_mapping)
        melted_df[self.choice2_colname] = melted_df[self.choice2_colname].map(choice2_mapping)


        # merge predicted probabilities with original df
        merged_df = pd.merge(new_df, melted_df,
                             on=[self.group_colname, self.choice_colname, self.choice2_colname],
                             how='inner')

        # fill in out of domain values
        common_cols = list(set(original_new_df) & set(merged_df))
        final_df = pd.merge(merged_df, original_new_df, how='outer', on=common_cols)

        # Fill missing values with a specified value
        final_df.fillna(0, inplace=True)
        return(final_df)

    def compute_loss(self, new_df):
        """
        Compute the loss for the given data.

        Returns:
            float: The computed loss.
        """
        choice1_variable_matrix, choice2_variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(new_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        loss = self.model.loss_fn(self.model.coefs,
                                  choice1_variable_matrix,
                                  choice2_variable_matrix,
                                  counts_matrix,
                                  mask_matrix)
        return(float(loss))


class ConditionalLogisticRegressionEvaluator(DataTransformer):
    def __init__(self, model_path, params, training_df = None, validation_df = None):
        self.model = self.load_model(model_path)
        self.model.training_df = training_df
        self.validation_df = validation_df
        if not isinstance(self.model, ConditionalLogisticRegressor):
            raise TypeError("'model' must be a ConditionalLogisticRegressor object")
        super().__init__(self.model.training_df, self.model.input_variable_colnames, self.model.count_colname, self.model.input_group_colname, self.model.repeat_obs_colname, self.model.input_choice_colname, params)
        self.log_loss = None
        self.expected_log_loss = None

    def load_model(self, file_path):
        with open(file_path, 'rb') as file:
            model = dill.load(file)
        assert model.coefs is not None, "model is not trained"
        return(model)

    def calculate_log_loss(self):
        assert self.model.training_df is not None, 'No input training dataframe provided'
        variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(self.model.training_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        # Compute the loss on the training data
        loss = self.model.loss_fn(self.model.coefs,
                                  variable_matrix,
                                  counts_matrix,
                                  mask_matrix)
        return(loss)

    def calculate_expected_log_loss(self, fold_count=20):
        assert self.model.training_df is not None, 'No input training dataframe provided'
        variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(self.model.training_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        # Compute the expected loss on the training data
        self.model.variable_matrix = variable_matrix
        self.model.counts_matrix = counts_matrix
        self.model.nonrepeat_grp_matrix = nonrepeat_grp_matrix

        expected = self.model.cv_loss(fold_count, self.model.l2reg)
        e_loss = sum((1/fold_count) * np.array(expected))
        return(e_loss)

    def calculate_validation_log_loss(self):
        assert self.validation_df is not None, 'No input validation dataframe provided'
        variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(self.validation_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        loss = self.model.loss_fn(self.model.coefs,
                                  variable_matrix,
                                  counts_matrix,
                                  mask_matrix)
        return(loss)

    def compile_evaluation_results_df(self, calculate_validation_loss = False, calculate_expected_loss=False):
        result = {'training_annotation_type':[self.params.annotation_type],
                  'productivity':[self.params.productivity],
                  'motif_length_5_end':[self.params.left_nuc_motif_count],
                  'motif_length_3_end':[self.params.right_nuc_motif_count],
                  'motif_type':[self.params.motif_type],
                  'gene_weight_type':[self.params.gene_weight_type],
                  'upper_bound':[self.params.upper_trim_bound],
                  'lower_bound':[self.params.lower_trim_bound],
                  'insertion_bound':[self.params.insertions],
                  'model_type':[self.params.model_type],
                  'base_count_5end_length':[10],
                  'model_parameter_count':[len(self.model.coefs)]}

        results_df = pd.DataFrame(result)

        final = pd.DataFrame()

        if calculate_expected_loss is True:
            self.expected_log_loss = self.calculate_expected_log_loss()

            # add expected loss result
            e = results_df.copy()
            e['loss_type'] = 'Expected log loss across training data'
            e['log_loss'] = self.expected_log_loss

            final = pd.concat([final, e], axis = 0)

        if calculate_validation_loss is True:
            val_loss = self.calculate_validation_log_loss()
            val = results_df.copy()
            val['loss type'] = 'Log loss on validation data'
            val['log_loss'] = val_loss

            final = pd.concat([final, val], axis = 0)

        return(final)

class TwoStepConditionalLogisticRegressionEvaluator(TwoStepDataTransformer):
    def __init__(self, model_path, params, training_df = None, validation_df = None):
        self.model = self.load_model(model_path)
        self.model.training_df = training_df
        self.validation_df = validation_df
        if not isinstance(self.model, TwoStepConditionalLogisticRegressor):
            raise TypeError("'model' must be a TwoStepConditionalLogisticRegressor object")
        super().__init__(self.model.training_df, self.model.input_variable_colnames, self.model.input_choice1_variable_colnames, self.model.input_choice2_variable_colnames, self.model.count_colname, self.model.input_group_colname, self.model.repeat_obs_colname, self.model.input_choice_colname, self.model.input_choice2_colname, params)
        self.log_loss = None
        self.expected_log_loss = None

    def load_model(self, file_path):
        with open(file_path, 'rb') as file:
            model = dill.load(file)
        assert model.coefs is not None, "model is not trained"
        return(model)

    def calculate_log_loss(self):
        assert self.model.training_df is not None, 'No input training dataframe provided'
        choice1_variable_matrix, choice2_variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(self.model.training_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        # Compute the loss on the training data
        loss = self.model.loss_fn(self.model.coefs,
                                  choice1_variable_matrix,
                                  choice2_variable_matrix,
                                  counts_matrix,
                                  mask_matrix)
        return(loss)

    def calculate_expected_log_loss(self, fold_count=20):
        assert self.model.training_df is not None, 'No input training dataframe provided'
        choice1_variable_matrix, choice2_variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(self.model.training_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        # Compute the expected loss on the training data
        self.model.choice1_variable_matrix = choice1_variable_matrix
        self.model.choice2_variable_matrix = choice2_variable_matrix
        self.model.counts_matrix = counts_matrix
        self.model.nonrepeat_grp_matrix = nonrepeat_grp_matrix

        expected = self.model.cv_loss(fold_count, self.model.l2reg)
        e_loss = sum((1/fold_count) * np.array(expected))
        return(e_loss)

    def calculate_validation_log_loss(self):
        assert self.validation_df is not None, 'No input validation dataframe provided'
        choice1_variable_matrix, choice2_variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(self.validation_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        loss = self.model.loss_fn(self.model.coefs,
                                  choice1_variable_matrix,
                                  choice2_variable_matrix,
                                  counts_matrix,
                                  mask_matrix)
        return(loss)

    def compile_evaluation_results_df(self, calculate_validation_loss = False, calculate_expected_loss=False):
        result = {'training_annotation_type':[self.params.annotation_type],
                  'productivity':[self.params.productivity],
                  'motif_length_5_end':[self.params.left_nuc_motif_count],
                  'motif_length_3_end':[self.params.right_nuc_motif_count],
                  'motif_type':[self.params.motif_type],
                  'gene_weight_type':[self.params.gene_weight_type],
                  'upper_bound':[self.params.upper_trim_bound],
                  'lower_bound':[self.params.lower_trim_bound],
                  'insertion_bound':[self.params.insertions],
                  'model_type':[self.params.model_type],
                  'base_count_5end_length':[10],
                  'model_parameter_count':[len(self.model.coefs)]}

        results_df = pd.DataFrame(result)

        final = pd.DataFrame()

        if calculate_expected_loss is True:
            self.expected_log_loss = self.calculate_expected_log_loss()

            # add expected loss result
            e = results_df.copy()
            e['loss_type'] = 'Expected log loss across training data'
            e['log_loss'] = self.expected_log_loss

            final = pd.concat([final, e], axis = 0)

        if calculate_validation_loss is True:
            val_loss = self.calculate_validation_log_loss()
            val = results_df.copy()
            val['loss type'] = 'Log loss on validation data'
            val['log_loss'] = val_loss

            final = pd.concat([final, val], axis = 0)

        return(final)
