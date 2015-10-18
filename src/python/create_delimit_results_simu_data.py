#! /usr/bin/env python
import os
import commands

def run_delimit_on_data(input_tree_file, output_delimit_single_minbr_0_file, output_delimit_multi_minbr_0_file, output_delimit_single_minbr_default_file, output_delimit_multi_minbr_default_file):
	try:
		open(input_tree_file)
		
		if not os.path.exists(os.path.dirname(output_delimit_single_minbr_0_file)):
	    			os.makedirs(os.path.dirname(output_delimit_single_minbr_0_file))
		if not os.path.exists(os.path.dirname(output_delimit_multi_minbr_0_file)):
	    			os.makedirs(os.path.dirname(output_delimit_multi_minbr_0_file))
		if not os.path.exists(os.path.dirname(output_delimit_single_minbr_default_file)):
	    			os.makedirs(os.path.dirname(output_delimit_single_minbr_default_file))
		if not os.path.exists(os.path.dirname(output_delimit_multi_minbr_default_file)):
	    			os.makedirs(os.path.dirname(output_delimit_multi_minbr_default_file))

		delimit_single_minbr_0_call = "./delimit --ml_single --min_br 0 --tree_file " + input_tree_file + " --output_file foo"
		delimit_multi_minbr_0_call = "./delimit --ml_multi --min_br 0 --tree_file " + input_tree_file + " --output_file foo"
		delimit_single_minbr_default_call = "./delimit --ml_single --tree_file " + input_tree_file + " --output_file foo"
		delimit_multi_minbr_default_call = "./delimit --ml_multi --tree_file " + input_tree_file + " --output_file foo"

		(stat_single_minbr_0, output_single_minbr_0) = commands.getstatusoutput(delimit_single_minbr_0_call)
		(stat_multi_minbr_0, output_multi_minbr_0) = commands.getstatusoutput(delimit_multi_minbr_0_call)
		(stat_single_minbr_default, output_single_minbr_default) = commands.getstatusoutput(delimit_single_minbr_default_call)
		(stat_multi_minbr_default, output_multi_minbr_default) = commands.getstatusoutput(delimit_multi_minbr_default_call)

		delimit_single_minbr_0_out = open(output_delimit_single_minbr_0_file, 'w')
		delimit_multi_minbr_0_out = open(output_delimit_multi_minbr_0_file, 'w')
		delimit_single_minbr_default_out = open(output_delimit_single_minbr_default_file, 'w')
		delimit_multi_minbr_default_out = open(output_delimit_multi_minbr_default_file, 'w')

		delimit_single_minbr_0_out.write(output_single_minbr_0)
		delimit_multi_minbr_0_out.write(output_multi_minbr_0)
		delimit_single_minbr_default_out.write(output_single_minbr_default)
		delimit_multi_minbr_default_out.write(output_multi_minbr_default)

		delimit_single_minbr_0_out.close()
		delimit_multi_minbr_0_out.close()
		delimit_single_minbr_default_out.close()
		delimit_multi_minbr_default_out.close()
	except IOError:
		print "File not found: " + input_tree_file

set_names = ["Ne1e+05", "Ne1e+06", "Ne5e+05", "Ne10000"]

for set_name in set_names:
	for i in range(1,101):
		input_tree_file = "SimulB_C_trees/set_" + set_name + "/SimulB_C_tree_set_" + set_name + "." + str(i) + ".txt"
		output_delimit_single_minbr_0_file = "SimulB_C_delimit_single_minbr_0/set_" + set_name + "/delimit_results_set_" + set_name + "." + str(i) + ".txt"
		output_delimit_multi_minbr_0_file = "SimulB_C_delimit_multi_minbr_0/set_" + set_name + "/delimit_results_set_" + set_name + "." + str(i) + ".txt"
		output_delimit_single_minbr_default_file = "SimulB_C_delimit_single_minbr_default/set_" + set_name + "/delimit_results_set_" + set_name + "." + str(i) + ".txt"
		output_delimit_multi_minbr_default_file = "SimulB_C_delimit_multi_minbr_default/set_" + set_name + "/delimit_results_set_" + set_name + "." + str(i) + ".txt"
		run_delimit_on_data(input_tree_file, output_delimit_single_minbr_0_file, output_delimit_multi_minbr_0_file, output_delimit_single_minbr_default_file, output_delimit_multi_minbr_default_file)
