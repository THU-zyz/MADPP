import sys

import pandas as pd

sys.path.append('./script')

import os
import argparse
import utils
import msd_analysis as ma

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--post_processed_file", default='', help="the path to input file, which is the result '.xlsx' file of post-processing")
parser.add_argument("-r", "--ref_table_file", default='', help="the path to ref file")
parser.add_argument("-o", "--output_file", type=str, default='./', help="the path to output file")
parser.add_argument("-f", "--FC", type=str, default=1.25, help="set a significant threshold for fold change")
parser.add_argument("-s", "--significance_test_method", type=str, default='TTEST', help="the significance testing method")
parser.add_argument("-fdr", "--FDR", type=bool, default=True, help="whether to use FDR")
parser.add_argument("-l", "--labels", type=list, help="list of testing subgroups")
parser.add_argument("-c", "--config", type=str, default='', help="the path to config file, the program will give priority ")

args = parser.parse_args()

configs_dict = {}
configs_dict['post_processed_file'] = args.post_processed_file
configs_dict['output_file_path'] = args.output_file
configs_dict['ref_table_file'] = args.ref_table_file
configs_dict['FC'] = args.FC
configs_dict['significance_test_method'] = args.significance_test_method
configs_dict['FDR'] = args.FDR
configs_dict['labels'] = args.labels
configs_dict['config'] = args.config
'''
if  len(args.output_name) < 0:
    configs_dict['output_name'] = os.path.basename(args.input_file)
else:
    configs_dict['output_name'] = args.output_name
'''
configs_dict = utils.da_resetting_config(args.config, configs_dict)
utils.da_check_config(configs_dict)

split_datas = ma.split_result_data(configs_dict['post_processed_file'],
                                   configs_dict['labels'],
                                   configs_dict['output_file_path'])

QC_df = pd.read_excel(configs_dict['post_processed_file'], sheet_name='QCs_TIC', index_col=0)
ma.TIC_map(QC_df, 'QCs', os.path.join(configs_dict['output_file_path'],'result'))

ma.MSDataAnalysis(input_path = configs_dict['output_file_path'],
                  output_path = os.path.join(configs_dict['output_file_path'],'result'),
                  significance_method = configs_dict['significance_test_method'],
                  isFDRcorrection = configs_dict['FDR'],
                  FC = configs_dict['FC'],
                  ref_table_file = configs_dict['ref_table_file'])

