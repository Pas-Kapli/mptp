#! /usr/bin/env python

import commands
import time

def evaluate(treeFile, rooted):
	cmd_multi = './delimit --ptp_multi --tree_file ' + treeFile + ' --output_file foo'
	cmd_single = './delimit --ptp_single --tree_file ' + treeFile + ' --output_file foo'
	cmd_ptp_rooted = './PTP/PTP.py -t ' + treeFile + ' -p -minbr 0 -o output -pvalue 1'
	cmd_ptp_unrooted = './PTP/PTP.py -t ' + treeFile + ' -p -minbr 0 -o output -pvalue 1 -r'

	if (rooted):
		programs = [cmd_multi, cmd_single, cmd_ptp_rooted]
		cmd_ptp = cmd_ptp_rooted
	else:
		programs = [cmd_multi, cmd_single, cmd_ptp_unrooted]
		cmd_ptp = cmd_ptp_unrooted

	scores = {}
	times = {}

	print "Testing " + treeFile + "..."

	# cmd_ptp:
	ts = time.time()
	( stat, output ) = commands.getstatusoutput(cmd_ptp)
	te = time.time()
	times['ptp'] = te-ts
	#print output
	left = output.find("MAX logl: ")
	right = output[left+10:].find("\n")
	score = output[left+10:right+left+10]
	scores['ptp'] = score

	# cmd_multi:
	ts = time.time()
	( stat, output ) = commands.getstatusoutput(cmd_multi)
	te = time.time()
	times['multi'] = te-ts
	#print output
	left = output.find("Best score found single: ")
	right = output[left+25:].find("\n")
	score = output[left+25:right+left+25]
	scores['multi'] = score

	# cmd_single:
	ts = time.time()
	( stat, output ) = commands.getstatusoutput(cmd_single)
	te = time.time()
	times['single'] = te-ts
	#print output
	left = output.find("Best score found single: ")
	right = output[left+25:].find("\n")
	score = output[left+25:right+left+25]
	scores['single'] = score
	
	print 'scores: '
	print scores
	print 'times: '
	print times
	print '\n'

	return scores

def compare_rooted():
	with open('tree_names_rooted') as f_rooted:
		content = f_rooted.read().splitlines()
	#gnuplotOut = open('workfile', 'w')
	for i in range (0, len(content)):
		scores = evaluate('trees/' + content[i], True)
		#gnuplotOut.write(str(i) + ' ' + scores['ptp'] + ' ' + scores['multi'] + ' ' + scores['single'] + '\n')
		#print evaluate('trees/' + name)
	#gnuplotOut.close()
	#commands.getstatusoutput('gnuplot plotscript')
	f_rooted.close()

def compare_unrooted():
	with open('tree_names_unrooted') as f_unrooted:
		content = f_unrooted.read().splitlines()
	#gnuplotOut = open('workfile', 'w')
	for i in range (0, len(content)):
		scores = evaluate('trees/' + content[i], False)
		#gnuplotOut.write(str(i) + ' ' + scores['ptp'] + ' ' + scores['multi'] + ' ' + scores['single'] + '\n')
		#print evaluate('trees/' + name)
	#gnuplotOut.close()
	#commands.getstatusoutput('gnuplot plotscript')
	f_unrooted.close()

compare_unrooted()
compare_rooted()
