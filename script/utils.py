# -*- coding: utf-8 -*-
import os
import sys
import pandas
def string_to_bool(value):
    if value.lower() == 'false':
        return False
    elif value.lower() == 'true':
        return True
    else:
        raise ValueError("Invalid literal for boolean: '{}'".format(value))
    
def dl_resetting_config(config_path, configs_dict):
    if len(config_path)<1:
        return configs_dict
    if not os.path.exists(config_path):
        sys.exit(f"Config file not found in {config_path}, please check the config file path")
    configs_file = open(config_path,'r')
    for line in configs_file.readlines():
        if len(line.split(' '))>2:
            configs_dict[line.split(' ')[0]] = [i.strip('\n') for i in line.split(' ')[1:]]
        else:
            configs_dict[line.split(' ')[0]] = line.split(' ')[1].strip('\n')
    configs_dict['sample_blank_ratio'] = float(configs_dict['sample_blank_ratio'])
    configs_dict['sample_blank_ratio_passing_rate'] = float(configs_dict['sample_blank_ratio_passing_rate'])
    configs_dict['qc_rsd'] = float(configs_dict['qc_rsd'])
    configs_dict['qc_rsd_passing_rate'] = float(configs_dict['qc_rsd_passing_rate'])
    if type(configs_dict['blacklist']) == type('str'):
        configs_dict['blacklist'] = [configs_dict['blacklist']]
    configs_dict['sample_blank_filter']=string_to_bool(configs_dict['sample_blank_filter'])
    configs_dict['qc_rsd_filter'] = string_to_bool(configs_dict['qc_rsd_filter'])
    return configs_dict

def dl_check_config(configs_dict):
    qc_rsd_range = [0,100]
    sample_blank_ratio_range=[0,20]
    replace_na_options = ['1k', 'half_min']
    neccessary_config = ['input_file_path']
    for config in configs_dict.keys():
        if len(str(configs_dict[config]))<1 and config in neccessary_config:
            print(f'Please enter nesccessary correct {config}')
            sys.exit()

    if configs_dict['replace_na_method']  not in replace_na_options:
        if configs_dict['replace_na_method'].isdigit():
            pass
        else:
            sys.exit(f'Please enter correct replace_na value, options: {replace_na_options} or other numbers')

    if configs_dict['qc_rsd']< qc_rsd_range[0] or configs_dict['qc_rsd'] > qc_rsd_range[1]:
        sys.exit(f'Input qc_rsd is out range of {qc_rsd_range}')

    if configs_dict['sample_blank_ratio']< sample_blank_ratio_range[0] or \
            configs_dict['sample_blank_ratio'] > sample_blank_ratio_range[1]:
        sys.exit(f'Input sample_blank_ratio is out range of {sample_blank_ratio_range}')

    print('Config loaded')

def da_resetting_config(config_path, configs_dict):
    if len(config_path)<1:
        return configs_dict
    if not os.path.exists(config_path):
        exit(f"Config file not found in {config_path}, please check the config file path")
    configs_file = open(config_path,'r')
    for line in configs_file.readlines():
        if len(line.split(' '))>2:
            configs_dict[line.split(' ')[0]] = [i.strip('\n') for i in line.split(' ')[1:]]
        else:
            configs_dict[line.split(' ')[0]] = line.split(' ')[1].strip('\n')
    configs_dict['FC'] = float(configs_dict['FC'])

    return configs_dict

def da_check_config(configs_dict):
    for config in configs_dict.keys():
        if len(str(configs_dict[config]))<1:
            print(f'Please enter correct {config}')
            sys.exit()
    significance_test_options = ['TTEST', 'MWU']
    if configs_dict['significance_test_method']  not in significance_test_options:
        sys.exit(f'Please enter correct replace_na value, options: {significance_test_options}')
    if len(configs_dict['labels']) < 2:
        sys.exit('Please enter correct labels of subgroups')
    print('Config loaded')