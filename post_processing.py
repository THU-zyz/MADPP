# -*- coding: utf-8 -*-
import sys
sys.path.append('./script')

import os
import argparse
import utils
import data_loader

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", default='', help="the path to input file, which contains raw batch files,"
                                                           "injection information and concentration reference table")
parser.add_argument("-o", "--output_file", type=str, default='./', help="the path to output file")
parser.add_argument("-r", "--replace_na_method", type=str, default='1k', help="the method of NA replacing, options: '1k', 'half_min'")
parser.add_argument("-s", "--sample_blank_ratio", type=float, default=1.5, help="the ratio of sample to blank")
parser.add_argument("-sp", "--sample_blank_ratio_passing_rate", type=float, default=0.8, help="proportion of batches with normal sb value")
parser.add_argument("-sf", "--sample_blank_filter", type=bool, default=False, help="If True, metabolites with unqualified s:b ratio will be removed")
parser.add_argument("-b", "--blacklist", type=str, default='', help="list of unwanted metabolites")
parser.add_argument("-q", "--qc_rsd", type=float, default=20, help="the specified rsd of the qc")
parser.add_argument("-qp", "--qc_rsd_passing_rate", type=float, default=0.8, help="proportion of batches with normal rsd")
parser.add_argument("-qf", "--qc_rsd_filter", type=bool, default=False, help="If True, metabolites with unqualified qc rsd will be removed")


#parser.add_argument("-n", "--output_name", type=str, default='', help="the name of output file, default value is name of input file")
parser.add_argument("-c", "--config", type=str, default='', help="the path to config file, the program will give priority ")

args = parser.parse_args()
configs_dict = {}
configs_dict['input_file_path'] = args.input_file
configs_dict['output_file_path'] = args.output_file
configs_dict['replace_na_method'] = args.replace_na_method
configs_dict['blacklist'] = args.blacklist
configs_dict['sample_blank_ratio'] = args.sample_blank_ratio
configs_dict['sample_blank_ratio_passing_rate'] = args.sample_blank_ratio_passing_rate
configs_dict['sample_blank_filter'] = args.sample_blank_filter

configs_dict['qc_rsd'] = args.qc_rsd
configs_dict['qc_rsd_passing_rate'] = args.qc_rsd_passing_rate
configs_dict['qc_rsd_filter'] = args.qc_rsd_filter


'''
if  len(args.output_name) < 0:
    configs_dict['output_name'] = os.path.basename(args.input_file)
else:
    configs_dict['output_name'] = args.output_name
'''

configs_dict = utils.dl_resetting_config(args.config, configs_dict)
utils.dl_check_config(configs_dict)

print(f"The config settings are listed below: ")
for key, value in configs_dict.items():
    print(f'{key}: {value}')
_check = input("Please confirm whether to continue[y/n]: ")
if not ((_check == "y") | (_check == "Y")):
    exit("Program terminated.")

#Loader = data_loader.Loader(configs_dict)
Processor = data_loader.Processor(configs_dict)

