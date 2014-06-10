import os, string, sys, re, subprocess, time

def wait(run_names):
	if not hasattr(run_names, '__iter__'):
		run_names = [run_names]
	while True:
		jlist = subprocess.Popen('jlist', shell=True, stdout=subprocess.PIPE).communicate()[0]
		if not any([(' '+run_name+' ' in jlist) for run_name in run_names]): break
		time.sleep(1)

