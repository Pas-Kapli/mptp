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

set_names = ["Ne10000", "Ne100000", "Ne500000", "Ne1000000"]

for set_name in set_names:
	for i in range(1,101):
		input_tree_file = "similar_to_GMYC/15-08-2015.16-40/set_BIRTH0.27_" + set_name + "/rooted.RAxML_result.inferred.simulated_set_BIRTH0.27_" + set_name + "_" + str(i) + ".phy"
		output_delimit_single_minbr_0_file = "similar_to_GMYC_delimit_single_minbr_0/set_" + set_name + "/delimit_results_set_" + set_name + "." + str(i) + ".txt"
		output_delimit_multi_minbr_0_file = "similar_to_GMYC_delimit_multi_minbr_0/set_" + set_name + "/delimit_results_set_" + set_name + "." + str(i) + ".txt"
		output_delimit_single_minbr_default_file = "similar_to_GMYC_delimit_single_minbr_default/set_" + set_name + "/delimit_results_set_" + set_name + "." + str(i) + ".txt"
		output_delimit_multi_minbr_default_file = "similar_to_GMYC_delimit_multi_minbr_default/set_" + set_name + "/delimit_results_set_" + set_name + "." + str(i) + ".txt"
		run_delimit_on_data(input_tree_file, output_delimit_single_minbr_0_file, output_delimit_multi_minbr_0_file, output_delimit_single_minbr_default_file, output_delimit_multi_minbr_default_file)
