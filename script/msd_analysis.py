# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import statistics
import math
from scipy import stats
import copy
from matplotlib import pyplot as plt
import seaborn as sns
def TIC_map(data_TIC, name, output_path):
    _h = 0.20 * data_TIC.shape[1]
    ticplot = plt.figure(figsize=(24, max(10, _h)))
    _left = data_TIC.shape[1] * [0]
    for i in range(0, data_TIC.shape[0]):
        plt.barh(data_TIC.columns, data_TIC.iloc[i, :], left=_left)
        _left = _left + data_TIC.iloc[i, :]
    plt.xlabel('TIC ratio')
    plt.ylabel("sample")
    plt.savefig(os.path.join(output_path,f"{name}_TIC.png"))
    plt.close(plt.get_fignums()[-1])

def split_result_data(post_processed_file, labels, output_path):
    data = pd.read_excel(post_processed_file,sheet_name='TIC',index_col = 0)
    split_datas = {}
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    split_res_path = os.path.join(output_path, 'subgroup_data')
    fig_res_path = os.path.join(output_path, 'result')
    if not os.path.exists(split_res_path):
        os.mkdir(split_res_path)
    if not os.path.exists(fig_res_path):
        os.mkdir(fig_res_path)
    for label in labels:
        data_TIC = data[data.columns[data.columns.str.contains(label)]].copy()
        exwriter = pd.ExcelWriter(os.path.join(split_res_path , f'{label}.xlsx'))
        data_TIC.to_excel(exwriter)
        split_datas[label]=data_TIC
        exwriter.close()
        TIC_map(data_TIC, label, fig_res_path)
    return split_datas



def get_annotations_from_ref(ref_path):
    ref = pd.DataFrame(pd.read_excel(ref_path, index_col=0))
    annotation_df = pd.DataFrame()
    annotation_df.index = ref.index
    if "HMDB" in ref.columns:
        annotation_df.insert(annotation_df.shape[1], column="HMDB", value=ref["HMDB"])
    else:
        annotation_df.insert(annotation_df.shape[1], column="HMDB", value="")
    if "PubChem" in ref.columns:
        annotation_df.insert(annotation_df.shape[1], column="PubChem", value=ref["PubChem"])
    else:
        annotation_df.insert(annotation_df.shape[1], column="PubChem", value="")
    if "KEGG" in ref.columns:
        annotation_df.insert(annotation_df.shape[1], column="KEGG", value=ref["KEGG"])
    else:
        annotation_df.insert(annotation_df.shape[1], column="KEGG", value="")
    if "Metaboanalyst_Hit" in ref.columns:
        annotation_df.insert(annotation_df.shape[1], column="Metaboanalyst_Hit", value=ref["Metaboanalyst_Hit"])
    else:
        annotation_df.insert(annotation_df.shape[1], column="Metaboanalyst_Hit", value="")

    return annotation_df

def draw_volcano(statistics_result_df, FC, pairs, _FDRstr, output_path):
    log2FC_threshold = math.log2(FC)
    # 火山图 # plsda图
    for key2 in pairs:
        plt.figure(figsize=(15, 18))  # 定义整张图的大小
        plt.axvline(x=log2FC_threshold, ls=":", color="grey", zorder=0)
        plt.axvline(x=-log2FC_threshold, ls=":", color="grey", zorder=0)
        plt.axhline(y=-math.log10(0.05), ls=":", color="grey", zorder=0)

        _x = "log2(FC):" + key2
        _y = "-log10(" + _FDRstr + "):" + key2
        _change = "Change(" + _FDRstr + "):" + key2

        ax = sns.scatterplot(x=_x, y=_y,
                             hue=_change,
                             hue_order=('UP', 'NORMAL', 'DOWN'),
                             palette=("#E41A1C", "grey", "#377EB8"),
                             data=statistics_result_df,
                             alpha=0.7,
                             legend=False
                             )

        for k in range(0, statistics_result_df.shape[0]):
            _index = statistics_result_df.index[k]
            if statistics_result_df.loc[_index, _change] != "NORMAL":
                _posx = statistics_result_df.loc[_index, _x]
                _posy = statistics_result_df.loc[_index, _y]
                if statistics_result_df.loc[_index, _change] == "UP":
                    plt.annotate(statistics_result_df.index[k], xy=(_posx, _posy), xytext=(_posx + 0.05, _posy + 0.05),
                                    color=(0.894, 0.102, 0.110), fontsize=8)
                elif statistics_result_df.loc[_index, _change] == "DOWN":
                    plt.annotate(statistics_result_df.index[k], xy=(_posx, _posy), xytext=(_posx + 0.05, _posy + 0.05),
                                    color=(0.2157, 0.494, 0.7216), fontsize=8)

        ax.set_ylabel("-log10(" + _FDRstr + ")", fontweight='bold')
        ax.set_xlabel('log2(FC)', fontweight='bold')
        ax.set_title("Volcano plot:" + key2)
        plt.savefig(os.path.join(output_path,f"{key2.replace('|','-')}_volcano.png"))
        plt.close(plt.get_fignums()[-1])
    print("Finished.")

def MSDataAnalysis(input_path, output_path, significance_method, isFDRcorrection, FC, ref_table_file):

    input_datas={}
    labels=[]
    pairs_dir = os.path.join(input_path, 'subgroup_data')
    for file_name in os.listdir(pairs_dir):
        data_path = os.path.join(pairs_dir, file_name)
        input_datas[file_name.replace('.xlsx','')] = pd.read_excel(os.path.join(data_path),index_col=0)
        labels.append(file_name.replace('.xlsx',''))

    if len(labels) >= 2:
        combination = int(
            math.factorial(len(labels)) / math.factorial(2) / math.factorial(len(labels) - 2))  # 计算组合数

    subgroup_mean = {}
    for key in input_datas.keys():
        subgroup_mean[key] = input_datas[key].mean(axis=1)

    subgroup_FC = {}
    subgroup_log2FC = {}
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            subgroup_FC[labels[i] + "|" + labels[j]] = subgroup_mean[labels[i]] / subgroup_mean[labels[j]]
    for key in subgroup_FC.keys():
        subgroup_log2FC[key] = np.log2(subgroup_FC[key])

    subgroup_pvalue = {}
    subgroup_log10pvalue = {}

    if significance_method == "TTEST":
        for i in range(0, len(labels)):
            for j in range(i + 1, len(labels)):
                _pvalue = []
                for _index in input_datas[labels[i]].index:
                    _pvalue.append(stats.ttest_ind(
                        input_datas[labels[i]].loc[_index, :],
                        input_datas[labels[j]].loc[_index, :],
                        alternative="two-sided", equal_var=False)[1])
                subgroup_pvalue[labels[i] + "|" + labels[j]]=_pvalue
    elif significance_method == "MWU":
        for i in range(0, len(labels)):
            for j in range(i + 1, len(labels)):
                _pvalue = []
                for _index in input_datas[labels[i]].index:
                    _pvalue.append(stats.mannwhitneyu(
                        input_datas[labels[i]].loc[_index, :],
                        input_datas[labels[j]].loc[_index, :],
                        alternative="two-sided")[1])
                subgroup_pvalue[labels[i] + "|" + labels[j]] = _pvalue
    else:
        exit("Invalid significant test arguments received.")
    for keys in subgroup_pvalue.keys():
        subgroup_log10pvalue[keys] = -np.log10(subgroup_pvalue[keys])
    labels2 = list(subgroup_pvalue.keys())

    # 计算FDR
    subgroup_FDR = copy.deepcopy(subgroup_pvalue)
    subgroup_log10FDR = {}
    for key2 in subgroup_FDR.keys():
        _sortindex = np.argsort(subgroup_FDR[key2])
        for j in range(0, len(subgroup_pvalue[key2])):
            subgroup_FDR[key2][_sortindex[j]] = subgroup_FDR[key2][_sortindex[j]] * len(subgroup_FDR[key2]) / (j + 1)
    for key2 in subgroup_FDR.keys():
        subgroup_log10FDR[key2] = -np.log10(subgroup_FDR[key2])

    # 判断是否用FDR还是p-value进行分类
    if isFDRcorrection:
        _standard = subgroup_FDR
    else:
        _standard = subgroup_pvalue
    subgroup_change = {}
    for key2 in labels2:
        subgroup_change[key2] = []
        for j in range(0, len(_standard[key2])):
            if _standard[key2][j] <= 0.05:
                if subgroup_FC[key2][j] >= FC:
                    subgroup_change[key2].append("UP")
                elif subgroup_FC[key2][j] <= 1/FC:
                    subgroup_change[key2].append("DOWN")
                else:
                    subgroup_change[key2].append("NORMAL")
            else:
                subgroup_change[key2].append("NORMAL")

    subgroup_mean_df = pd.DataFrame()

    for key in subgroup_mean.keys():
        subgroup_mean_df.insert(subgroup_mean_df.shape[1], column="Mean:" + key, value=subgroup_mean[key])

    sig_df = pd.DataFrame()
    if isFDRcorrection:
        _FDRstr = "FDR"
    else:
        _FDRstr = "p-value"

    for key2 in labels2:
        sig_df.insert(sig_df.shape[1], column="FC:" + key2, value=subgroup_FC[key2])
        sig_df.insert(sig_df.shape[1], column="FC:" + key2.split('|')[1]+'|'+key2.split('|')[0],
                      value=1/subgroup_FC[key2])

        sig_df.insert(sig_df.shape[1], column="log2(FC):" + key2, value=subgroup_log2FC[key2])
        sig_df.insert(sig_df.shape[1], column="log2(FC):" + key2.split('|')[1]+'|'+key2.split('|')[0],
                      value=-subgroup_log2FC[key2])

        sig_df.insert(sig_df.shape[1], column="p-value:" + key2, value=subgroup_pvalue[key2])
        sig_df.insert(sig_df.shape[1], column="-log10(p-value):" + key2, value=subgroup_log10pvalue[key2])
        sig_df.insert(sig_df.shape[1], column="FDR:" + key2, value=subgroup_FDR[key2])
        sig_df.insert(sig_df.shape[1], column="-log10(FDR):" + key2, value=subgroup_log10FDR[key2])
        sig_df.insert(sig_df.shape[1], column="Change(" + _FDRstr + "):" + key2, value=subgroup_change[key2])

    annotation_df = get_annotations_from_ref(ref_table_file)
    annotation_df = annotation_df.loc[subgroup_mean_df.index, :]
    statistics_result_df = pd.concat([annotation_df, subgroup_mean_df, sig_df], axis=1)

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    exwriter = pd.ExcelWriter(os.path.join(output_path, 'statistics_result.xlsx'))

    statistics_result_df.to_excel(exwriter, sheet_name='FC&p-value')

    for key2 in labels2:
        _change = "Change(" + _FDRstr + "):" + key2

        _subgroup_res_df = statistics_result_df[(statistics_result_df[_change] == "UP") | (statistics_result_df[_change] == "DOWN")].copy()

        _subgroup_res_df.to_excel(exwriter, sheet_name="result_" + key2)

    draw_volcano(statistics_result_df, FC, labels2, _FDRstr, output_path)



    exwriter.close()


