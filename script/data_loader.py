# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import pandas as pd

class Loader():
    def __init__(self, config_dict):
        self.config_dict = config_dict
        self.input_file_path = config_dict['input_file_path']
        self.blacklist = self.config_dict['blacklist']
        self.replace_na_method = self.config_dict['replace_na_method']
        self.sub_folder_names = ['concentration_table', 'injection_information', 'raw_batch']
        self.input_sub_folders = os.listdir(self.input_file_path)
        for sub_folder_name in self.sub_folder_names:
            if sub_folder_name not in self.input_sub_folders:
                exit(f'Missing or incorrect name of {sub_folder_name} in input file')
        self.conc_table_folder_path = os.path.join(self.input_file_path, 'concentration_table')
        self.injection_information_folder_path = os.path.join(self.input_file_path, 'injection_information')
        self.raw_batch_folder_path = os.path.join(self.input_file_path, 'raw_batch')

        self.conc_table = pd.DataFrame(pd.read_excel(os.path.join(self.conc_table_folder_path, 'ref.xlsx'), index_col=0))

        self.injection_sheet_check()
        self.get_dilution_set()
        self.combine_df()
        self.remove_metabolite()
        self.get_na_position(self.m_removed_combined)
        self.replace_NA()
        self.merge_concentration()


    def injection_sheet_check(self):
        '''
        Description: 检查injection_sheet的名称与 raw batch 文件是否对应。如果batch_name和sheet_name对应不上，则报错。
        -----------
        Parameter:
        ----------
        self.injection_information_folder_path
        self.raw_batch_folder_path
        return:
        -------
        self.batch_names
        '''
        injection_info_file_path = os.path.join(self.injection_information_folder_path, "injection_information.xlsx")
        # Retrieve the names of all worksheets from the specified Excel file.
        sheet_names = pd.ExcelFile(injection_info_file_path).sheet_names
        self.batch_names = [i.split('.csv')[0] for i in os.listdir(self.raw_batch_folder_path)]
        # 检查 self.batch_names 中的每个 batch_name 是否在 sheet_names 中
        for batch_name in self.batch_names:
            if batch_name not in sheet_names:
                exit(
                    f'Missing sheet name: {batch_name} in injection_information.xlsx, which is not matched with raw batch file name')
        # 检查 sheet_names 中的每个 sheet_name 是否在 self.batch_names 中
        for sheet_name in sheet_names:
            if sheet_name not in self.batch_names:
                exit(
                    f'Missing raw batch file for sheet name: {sheet_name} in injection_information.xlsx')

    def get_dilution_set(self):
        # 遍历conc_table的列名，筛选满足特定条件的列名：1.检查列名的小写是否为'ori'； 2.检查用strip去除列名中所有的x(不区分大小写)后，剩余字符串是否只包含数字。
        self.dilution_set = [i for i in self.conc_table.columns.values if i.lower() == 'ori' or str(i.strip('x').strip('X')).isdigit()]
        # dilution_set应该只有2个。
        if len(self.dilution_set)>2:
            exit("Too many dilution, please check the number of dilution in conc table, no more than 2 dilution")
        elif len(self.dilution_set)<2:
            exit("Dilution sets is less than 2, please check dilution in conc table")

    def ref_metabolites_check(self, batch_df, conc_table):
        batch_index = batch_df.index.values
        conc_table_index = conc_table.index.values
        not_in_table = [i for i in batch_index if i not in conc_table_index]
        not_in_batch = [i for i in conc_table_index if i not in batch_index]
        return not_in_table+not_in_batch

    def check_sample_names(self, sample_names, batch_combined_df):
        # 当某一个样本在batch中出现却未在sample_names中出现则记录到missing sample中
        # 当某一个样本在sample_names(进样顺序)中出现却未在batch中出现则记录到unfound sample中，这种情况将报错。
        self.missing_sample = [i for i in batch_combined_df.columns.values if i not in sample_names]
        print(f"Warning, these samples are not recored in injection_information.xlsx: {self.missing_sample.__str__()})")
        self.unfound_sample = [i for i in sample_names if i not in batch_combined_df.columns.values]
        if len(self.unfound_sample)>0:
            sys.exit(f'These samples in injection_information.xlsx could not be found in raw batch data: {self.unfound_sample.__str__()}')

    def combine_df(self):
        batch_files = os.listdir(self.raw_batch_folder_path)
        self.n_batch = len(batch_files)
        print("Batch file list: " + sorted(batch_files).__str__())
        # input用于获取用户的输入，即确认上述的信息是否正确，根据输入的值，决定后续数据处理过程。
        _check = input("Please confirm the information. Correct?[y/n]: ")
        if not ((_check == "y") | (_check == "Y")):
            exit("Program terminated.")

        raw_batch_combined = pd.DataFrame()
        sample_names = np.array([])

        self.dismatched_metabolites = {}
        for batch_file in batch_files:
            batch_name = batch_file.split('.csv')[0]
            tmp_batch_df = pd.read_csv(os.path.join(self.raw_batch_folder_path, batch_file),
                                       sep=",", index_col=0, encoding="UTF-8")
            # 去除含有Sample Type和Acquisition Date & Time的行
            tmp_batch_df = tmp_batch_df.loc[
                tmp_batch_df.index[~tmp_batch_df.index.str.contains("Sample Type|Acquisition Date & Time")],]
            # 把batch name加到每列列名之前，这样数据合并后，就知道每列数据来自的batch了：
            for j in range(0, tmp_batch_df.shape[1]):
                tmp_batch_df.columns.values[j] = "_".join([batch_name, tmp_batch_df.columns.values[j]])
            # check whether the metabolites in ref.xlsx is matched with those in batch_file

            self.dismatched_metabolites[batch_name] = self.ref_metabolites_check(tmp_batch_df, self.conc_table)
            if len(self.dismatched_metabolites[batch_name])>0:
                exit(f'Metabolites: {self.dismatched_metabolites[batch_name]} are not matched between {batch_file}.csv and ref.xlsx')

            raw_batch_combined = pd.concat([raw_batch_combined, tmp_batch_df], axis=1)

            # get sample names from injection information file
            tmp_sample_names = pd.read_excel(os.path.join(self.injection_information_folder_path,"injection_information.xlsx"),
                                             sheet_name=batch_name,header=None)
            for j in range(0, tmp_sample_names.shape[0]):
                tmp_sample_names.iloc[:, 0].values[j] = "_".join([batch_name, tmp_sample_names.iloc[:, 0].values[j]])
            sample_names = np.append(sample_names, tmp_sample_names.iloc[:, 0].values)


        self.check_sample_names(sample_names, raw_batch_combined)

        # remove irrelevant QC:
        # 筛选出所有包含子字符串QC的元素，讲这些元素转换为NumPy数组。
        QC_names = np.array([sample_name for sample_name in sample_names if "QC" in sample_name])
        QC_columns = raw_batch_combined[QC_names]
        raw_batch_combined = raw_batch_combined[
            raw_batch_combined.columns[~raw_batch_combined.columns.str.contains("QC")]]
        raw_batch_combined = pd.concat([raw_batch_combined, QC_columns], axis=1)
        raw_batch_combined = raw_batch_combined.astype(float)
        raw_batch_combined.index.name = "Metabolite"

        self.raw_batch_combined = raw_batch_combined
        self.sample_names = sample_names
        self.QC_names = QC_names

    def remove_unwanted(self, conc_table, raw_batch_combined):
        # 根据black_list中给出的代谢物名称从raw_batch_combined中去除相应代谢物。
        undetected_removed_df = raw_batch_combined.copy()
        for metabolite in self.blacklist:
            undetected_removed_df.drop(undetected_removed_df.index[undetected_removed_df.index.str.contains(metabolite)],
                                       axis=0, inplace=True)
            conc_table.drop(conc_table.index[conc_table.index.str.contains(metabolite)], axis=0, inplace=True)
            print("Remove metabolite: " + metabolite)
        return conc_table, undetected_removed_df

    def remove_undeteced(self, concentration_table_df, unwanted_removed_df):
        ref = concentration_table_df
        unwanted_removed_df = unwanted_removed_df.loc[ref.index.values] # remove the metabolites that didn't appear in ref

        self.undetected = []
        if self.dilution_set[0].lower() == 'ori':
            ori_col = ref[self.dilution_set[0]].copy()
            dil_col = ref[self.dilution_set[1]].copy()
        else:
            ori_col = ref[self.dilution_set[1]].copy()
            dil_col = ref[self.dilution_set[0]].copy()
        for index in ref.index.values:
            if str(ori_col.loc[index]).lower() == "x" and str(dil_col.loc[index]).lower() == "x":
                sys.exit(f'Wrong concentration selection in {index} in ref.xlsx')
            elif str(ori_col.loc[index]).lower() != "x" and str(dil_col.loc[index]).lower() != "x":
                ref.drop(index, inplace=True)
                self.undetected.append(index)
                if index in unwanted_removed_df.columns.values:
                    unwanted_removed_df.drop(index, inplace=True)
                else:
                    pass

        return ref, unwanted_removed_df

    def remove_metabolite(self):
        unwanted_removed_conc_table, unwanted_removed_df = self.remove_unwanted(self.conc_table, self.raw_batch_combined)
        self.m_removed_conc_table, self.m_removed_combined = self. remove_undeteced(unwanted_removed_conc_table,
                                                                                         unwanted_removed_df)

    def get_na_position(self, df):# 要随着rsd和sb一起算吗
        self.is_na_df = pd.isna(df)
        self.na_positions = [(row, col) for row in self.is_na_df.index
                             for col in self.is_na_df.columns if self.is_na_df.at[row, col]]

    def replace_na_with_halfmin(self, row):
        half_min = row.min()/2
        return row.fillna(half_min)

    def replace_NA(self):
        tmp_df = self.m_removed_combined
        if self.replace_na_method == '1k':
            replace_value = 1000.0
            tmp_df.fillna(value=replace_value, inplace=True)
        elif self.replace_na_method == 'half_min':
            tmp_df = tmp_df.apply(self.replace_na_with_halfmin, axis=1)
        elif self.replace_na_method.isdigit():
            replace_value = float(self.replace_na_method)
            tmp_df.fillna(value=replace_value, inplace=True)
        else:
            exit('Please enter correct replace_na_method')
        # 将tmp_df中所有等于0的元素替换为1000。
        tmp_df.mask(tmp_df == 0, other=1000, inplace=True)
        self.na_m_removed_combined = tmp_df

    def merge_concentration(self):
        ref = self.m_removed_conc_table
        conc_merged = pd.DataFrame()
        ori_df = self.na_m_removed_combined
        self.actual_sample=pd.DataFrame() # the actual selected dilution of each metabolite
        dil_index = {}
        # 根据 ref 中稀释度 dil 的值是否为 ‘x’ 或 ‘X’ 来获取该稀释度的索引。
        for dil in sorted(self.dilution_set):
            dil_index[dil] = ref[np.logical_or(ref[dil].values == 'x', ref[dil].values == 'X')].index.values
        # 根据 dil_index 中的索引，从 ori_df 中选择对应的行，并过滤出包含 key 的列。
        for key in dil_index.keys():
            specified_conc_subdf = ori_df.loc[dil_index[key]].filter(like=key, axis=1)

            tmp_conc_merged = pd.DataFrame()
            self.actual_sample = np.append(self.actual_sample, specified_conc_subdf.columns.values)
            '''
            for column_name in specified_conc_subdf.columns.values:
                tmp_col = specified_conc_subdf[column_name]
                tmp_col.columns = [column_name.replace(f'_{key}', '')]
                tmp_conc_merged = pd.concat([tmp_conc_merged, tmp_col],axis=1)
                #tmp_conc_merged[column_name.replace(f'_{key}', '')] = specified_conc_subdf[column_name]
            '''
            tmp_cols = specified_conc_subdf.copy()
            tmp_cols.columns = [i.replace(f'_{key}', '') for i in tmp_cols.columns.values]
            tmp_conc_merged = pd.concat([tmp_conc_merged, tmp_cols], axis=1)
            '''
            if key != 'Ori':
                tmp_conc_merged = tmp_conc_merged.mul(float(key.strip('x').strip('X')))
            '''
            if len(conc_merged)>0:
                # ValueError: The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()
                check_point=False
                if (len(conc_merged.columns) == len(tmp_conc_merged.columns)):
                    check_point=(conc_merged.columns == tmp_conc_merged.columns).all()

                if (check_point):
                    conc_merged = pd.concat([conc_merged, tmp_conc_merged], axis=0)
                else:
                    # Analyze the correspondence between samples at different dilution factors
                    ld = []
                    ld1 = [i for i in tmp_conc_merged.columns.values if i not in conc_merged.columns.values]
                    ld2 = [i for i in conc_merged.columns.values if i not in tmp_conc_merged.columns.values]
                    for i in ld1:
                        for j in conc_merged.columns.values:
                            if "_".join(i.split('_')[:-1]) in j and i.split('_')[-1].isdigit():
                                ld.append(i)
                                break
                    for i in ld2:
                        for j in tmp_conc_merged.columns.values:
                            if "_".join(i.split('_')[:-1]) in j and i.split('_')[-1].isdigit():
                                ld.append(i)
                                break
                    tmp_ld = ld1+ld2
                    ld = [i for i in tmp_ld if i not in ld]

                    if len(ld)>0:
                        print(f"Dilution({key}) setting is not matched in {ld},")
                        _check = input("Please confirm the information: delete the redundant samples ensuring program continuing[y/n]: ")
                        if not ((_check == "y") | (_check == "Y")):
                            exit("Program terminated.")
                        else:
                            print('Redundant samples are deleted')

                    conc_merged = pd.concat([conc_merged.drop(columns=ld), tmp_conc_merged.drop(columns=ld)], axis=0)
            else:
                conc_merged = pd.concat([conc_merged, tmp_conc_merged], axis=0)
        conc_merged.index.name = "Metabolite"
        self.conc_merged_df = conc_merged


class Processor(Loader):
    def __init__(self, config_dict):
        super().__init__(config_dict)
        self.output_file_path = self.config_dict['output_file_path']
        self.sample_blank_ratio = self.config_dict['sample_blank_ratio']
        self.sample_blank_ratio_passing_rate = self.config_dict['sample_blank_ratio_passing_rate']
        self.sample_blank_filter = self.config_dict['sample_blank_filter']
        self.qc_rsd = self.config_dict['qc_rsd']
        self.qc_rsd_passing_rate = self.config_dict['qc_rsd_passing_rate']
        self.qc_rsd_filter = self.config_dict['qc_rsd_filter']
        self.sb_unqualified = np.array([])
        self.rsd_unqualified = np.array([])

        self.sb_rsd_check()
        self.batch_normalize()
        self.TIC_normalize()
        self.dump_data()
        self.report()

    def blank_check(self, data_rebuild, Blank_tag='Blank'):
        for batch_name in self.batch_names:
            batch_subdf = data_rebuild.filter(like=batch_name)
            batch_subdf = batch_subdf[batch_subdf.columns[~batch_subdf.columns.str.contains('QC')]]
            blank_values = batch_subdf[batch_subdf.columns[batch_subdf.columns.str.contains(Blank_tag)]]
            if len(blank_values.columns.values) == 0:
                _check = input(
                    f"Warning: Blank sample doesnt exist in {batch_name}, continue? [y/n]: ")
                if not ((_check == "y") | (_check == "Y")):
                    exit("Program terminated.")
                continue

    def sb_check(self, data_rebuild, Blank_tag='Blank'):
        sb_check_df = pd.DataFrame()
        self.batch_sbratio_df = pd.DataFrame()
        _sb_ratio_dict = {}

        self.blank_check(data_rebuild, Blank_tag='Blank')

        for i in range(0, data_rebuild.shape[0]):
            tmp_row_df = data_rebuild.iloc[i:i + 1, :]
            batchs_sb_ratio = np.array([])
            for batch_name in self.batch_names:
                batch_subdf = tmp_row_df.filter(like=batch_name)
                batch_subdf = batch_subdf[batch_subdf.columns[~batch_subdf.columns.str.contains('QC')]]
                blank_values = batch_subdf[batch_subdf.columns[batch_subdf.columns.str.contains(Blank_tag)]]

                s_values = batch_subdf[batch_subdf.columns[~batch_subdf.columns.str.contains(Blank_tag)]]

                blank_values_mean = blank_values.mean(axis=1)
                s_values_mean = s_values.mean(axis=1)

                _sb_ratio = s_values_mean.values / blank_values_mean.values
                batchs_sb_ratio = np.append(batchs_sb_ratio, _sb_ratio)
                _sb_ratio_dict[batch_name] = _sb_ratio

            _index = tmp_row_df.index.values[0]
            _df = pd.DataFrame(_sb_ratio_dict, index = [_index])
            self.batch_sbratio_df = pd.concat([self.batch_sbratio_df,_df])

            percentage_above = np.sum(batchs_sb_ratio > self.sample_blank_ratio) / len(batchs_sb_ratio)
            
            if percentage_above > float(self.sample_blank_ratio_passing_rate):
                sb_check_df = pd.concat([sb_check_df, tmp_row_df])
            else:
                if self.sample_blank_filter:
                    self.sb_unqualified = np.append(self.sb_unqualified, _index)
                else:
                    sb_check_df = pd.concat([sb_check_df, tmp_row_df])
                    self.sb_unqualified = np.append(self.sb_unqualified, _index)

        sb_check_df[sb_check_df.columns[~sb_check_df.columns.str.contains(Blank_tag)]]  # drop every Blank column
        return sb_check_df

    def rsd_check(self, data_rebuild):
        rsd_check_df = pd.DataFrame()
        self.batch_rsd_df = pd.DataFrame()
        for i in range(0, data_rebuild.shape[0]):
            tmp_row_df = data_rebuild.iloc[i:i+1,:]
            batchs_rsd_ratio = np.array([])
            _rsd_ratio_dict = {}
            for batch_name in self.batch_names:

                batch_subdf = tmp_row_df.filter(like=batch_name)
                batch_QC_subdf = batch_subdf[batch_subdf.columns[batch_subdf.columns.str.contains('QC')]]
                _rsd_ratio = float(100*(np.std(batch_QC_subdf.values)/np.mean(batch_QC_subdf.values)))
                batchs_rsd_ratio = np.append(batchs_rsd_ratio, _rsd_ratio)
                _rsd_ratio_dict[batch_name] = _rsd_ratio

            _index = tmp_row_df.index.values[0]
            _df = pd.DataFrame(_rsd_ratio_dict, index = [_index])
            self.batch_rsd_df = pd.concat([self.batch_rsd_df,_df])

            percentage_above = np.sum(batchs_rsd_ratio < self.qc_rsd) / len(batchs_rsd_ratio)
            if percentage_above > float(self.qc_rsd_passing_rate):
                rsd_check_df = pd.concat([rsd_check_df, tmp_row_df])
            else:
                if self.qc_rsd_filter:
                    self.rsd_unqualified = np.append(self.rsd_unqualified, _index)
                else:
                    rsd_check_df = pd.concat([rsd_check_df, tmp_row_df])
                    self.rsd_unqualified = np.append(self.rsd_unqualified, _index)
        return rsd_check_df

    def sb_rsd_check(self):
        # 去除含有MPA，std的列

        tmp_data_rebuild = self.conc_merged_df.copy()
        _df = tmp_data_rebuild[tmp_data_rebuild.columns[~tmp_data_rebuild.columns.str.contains("MNA|Std|MPA_|std_high|std_low")]]

        data_rebuild_sb = self.sb_check(_df)

        data_rebuild_rsd = self.rsd_check(data_rebuild_sb)

        self.data_rebuild_df = data_rebuild_rsd

    def injection_pair_check(self, inj_seq):
        dil_1 = inj_seq.loc[inj_seq[0].str.contains(self.dilution_set[0])][0].values
        dil_2 = inj_seq.loc[inj_seq[0].str.contains(self.dilution_set[1])][0].values
        if len(dil_1) != len(dil_2):
            exit(f"Error occurred: dilution tag mismatch detected: {dil_1} {dil_2}")

    def batch_normalize(self):
        # Calculate the mean peak area of all the QC samples in all batches: df.QC.all_mean
        df_qc = self.data_rebuild_df.filter(like="QC")
        # df_qc.set_index(self.data_rebuild_df.index, inplace=True)


        tmp_col = df_qc.mean(axis=1).to_frame()
        tmp_col.columns = ["mean_QC"]
        df_qc = pd.concat([df_qc,tmp_col],axis=1)

        # QC normalization
        df_final_s = pd.DataFrame(index=self.data_rebuild_df.index)
        df_final_qc = pd.DataFrame(index=self.data_rebuild_df.index)

        # Read all sheet names from the Excel file (batch file names)
        injection_info_file_path = os.path.join(self.injection_information_folder_path, "injection_information.xlsx")
        # self.bacth_normalized_df = pd.concat([df_final_s, self.data_rebuild_df.filter(like="QC")], axis=1)
        for batch_name in self.batch_names:
            print("************************")
            print(f"Processing: {batch_name}")
            # Sample injection sequence of each batch
            inj_seq = pd.read_excel(injection_info_file_path, sheet_name=batch_name, header=None)
            self.injection_pair_check(inj_seq)
            inj_seq = inj_seq[inj_seq[0].str.contains(self.dilution_set[0])]  # anchor dilution 0
            # 对inj_seq重新编号
            inj_seq = inj_seq.reset_index(drop=True)
            for i in range(inj_seq.shape[0]):
                inj_seq.loc[i][0] = batch_name + '_' + inj_seq.loc[i][0].replace('_' + self.dilution_set[0], '')

            # inj_seq = inj_seq[inj_seq.iloc[:, 0].isin(self.actual_sample)]

            # QC index in each batch

            qc_index = inj_seq[inj_seq[0].str.contains("QC")].index.tolist()
            batch_qc_num = len(qc_index)

            for j in range(batch_qc_num):
                if j != 0 and qc_index[j] - qc_index[j - 1] > 1:
                    qc1 = f"{inj_seq.iloc[qc_index[j - 1]][0]}"
                    qc2 = f"{inj_seq.iloc[qc_index[j]][0]}"
                    tmp_group_qc = pd.concat([df_qc[qc1], df_qc[qc2]], axis=1)

                    # Recognize samples between each pairs of QC, and get the corresponding name: (batchname)_(sample_name)
                    group_j = inj_seq.iloc[qc_index[j - 1] + 1:qc_index[j]][0].values
                    # Get the corresponding column index of the samples in the combined data
                    # group_j_col = group_j.map(p.data_rebuild_df.columns).fillna(-1).astype(int)
                    group_j_col = self.data_rebuild_df.loc[:, self.data_rebuild_df.columns.isin(group_j)].columns.values

                    if len(group_j_col) != len(group_j):
                        # If there are -1 values in group_index, there is a name match failure
                        print(f"Name match failure detected: {group_j[group_j_col == -1].tolist()}")
                        continue

                    # Calculate the mean value of QC samples in each group
                    group_qc_mean = tmp_group_qc.mean(axis=1)

                    # Calculate the normalization factor of each group
                    group_qc_factor = df_qc["mean_QC"] / group_qc_mean
                    print("QC block:")
                    print(f"QC1: {qc1}")
                    print(f"QC2: {qc2}")
                    for k in group_j:
                        print(f"SAMPLE: {k}")

                    # Sample QC_mean Normalization
                    df_s_qcm_normalized = self.data_rebuild_df[group_j_col] * group_qc_factor.values[:, np.newaxis]
                    df_final_s = pd.concat([df_final_s, df_s_qcm_normalized], axis=1)

        self.bacth_normalized_df = df_final_s
        self.bacth_normalized_df.index.name = "Metabolite"
        self.QCs_df = df_qc
        #self.bacth_normalized_df = pd.concat([df_final_s, self.data_rebuild_df.filter(like="QC")], axis=1)

    def TIC_normalize(self):
        self.TIC_normalized_df = self.bacth_normalized_df.apply(lambda col: col / col.sum(), axis=0)
        self.QCs_df = self.QCs_df.apply(lambda col: col / col.sum(), axis=0)
        _count = 1
        form_warning = pd.DataFrame(columns=["Sample", "Metabolites", "WarningMessage", "value"])
        data_TIC = self.TIC_normalized_df.copy()
        for i in range(0, data_TIC.shape[0]):
            for j in range(0, data_TIC.shape[1]):
                if data_TIC.iloc[i, j] > 0.1:  # 如果TIC占比超过10%
                    form_warning.loc[_count] = [data_TIC.columns[j], data_TIC.index[i], "High ratio of total peak area",
                                                data_TIC.iloc[i, j]]
                    _count = _count + 1
                elif data_TIC.iloc[i, j] < 0.0001:  # 如果TIC占比低于0.01%
                    form_warning.loc[_count] = [data_TIC.columns[j], data_TIC.index[i], "Low ratio of total peak area",
                                                data_TIC.iloc[i, j]]

        #self.form_warning = form_warning.sort_values("value", ascending=False, inplace=True)  # 对挑选出来的条目根据TIC的数值大小进行降序排序，为了输出的时候整洁一些

        self.form_warning = form_warning


    def dump_data(self):
        if not os.path.exists(self.output_file_path):
            os.mkdir(self.output_file_path)
        '''
        self.na_m_removed_combined.to_csv(os.path.join(self.output_file_path, "raw_data_combined.csv"), encoding="UTF-8")
        print(f"The combined file has been saved to {self.output_file_path}/raw_data_combined.csv")
        self.conc_merged_df.to_csv(os.path.join(self.output_file_path, "conc_merged.csv"), encoding="UTF-8")
        print(f'The concentration integrated file has been saved to {self.output_file_path}/conc_merged.csv')
        '''
        # 用于将各步骤中的关键数据存储到result.xlsx这一个文件中。这里用到了ExcelWriter函数，并指定engine='openpyxl'。
        exwriter = pd.ExcelWriter(os.path.join(self.output_file_path, "result.xlsx"), engine='openpyxl')


        self.na_m_removed_combined.sort_index(axis=1).to_excel(exwriter, sheet_name="raw_data_combined" )
        self.bacth_normalized_df.sort_index(axis=1).to_excel(exwriter, sheet_name="QC_filtered")
        self.TIC_normalized_df.sort_index(axis=1).to_excel(exwriter, sheet_name="TIC")
        self.form_warning.to_excel(exwriter, sheet_name="Warning")
        self.batch_rsd_df.sort_index(axis=1).to_excel(exwriter, sheet_name="RSD")
        self.batch_sbratio_df.sort_index(axis=1).to_excel(exwriter, sheet_name="sample-blank")
        self.QCs_df.sort_index(axis=1).to_excel(exwriter, sheet_name="QCs_TIC")

        exwriter.close()

        print(f"Done. the normalized data file has been saved to {self.output_file_path}/result.xlsx")

        pass

    def report(self):
        output ="Metabolite type: {}\n" \
                "Undetected metabolites: {}\n"\
                "Unqualified 'sample:blank' setting metabolites: {}\n" \
                "Unqualified 'RSD' setting metabolites: {}\n"\
                "Irrelevant samples in raw batch files: {}\n"\
                .format(
                len(self.TIC_normalized_df),
                self.undetected.__str__(),
                self.sb_unqualified,
                self.rsd_unqualified,
                self.missing_sample)
        print(output)
        f = open(os.path.join(self.output_file_path, 'report.txt'),'w')
        f.write(output)
        f.close()