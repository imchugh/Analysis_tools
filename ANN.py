# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 14:26:14 2016

@author: imchugh
"""

from ffnet import ffnet, tmlgraph
import pdb

def main(inputs_arr, target_arr, configs_dict):

    # Do basic checks
    if not inputs_arr.shape[0] == target_arr.shape[0]:
        raise Exception('Input and target arrays must have same first ' \
                        'dimension... Returning!')

    n_input_nodes = configs_dict['node_arch_list'][0]
    if not n_input_nodes == inputs_arr.shape[1]:
        raise Exception('Specified input node architecture (n = %s) ' \
                        'incompatible with passed input arrays... Returning!'
                        %str(n_input_nodes))

    n_target_nodes = configs_dict['node_arch_list'][-1]
    if len(target_arr.shape) == 1:
        sec_dim = 1
    else: 
        sec_dim = target_arr.shape[1]
    if not n_target_nodes == sec_dim:
        raise Exception('Specified target node architecture (n = %s) ' \
                        'incompatible with passed input arrays... Returning!'
                        %str(n_target_nodes))        

    # Generate network and train
    conec = tmlgraph(configs_dict['node_arch_list'])
    net = ffnet(conec)
    net.train_tnc(inputs_arr, target_arr, 
                  maxfun = configs_dict['iterations'], messages = 1)

    # Do testing if requested
    if configs_dict['test?']:
        return net.test(inputs_arr, target_arr)
    else:
        return net.call(inputs_arr)
