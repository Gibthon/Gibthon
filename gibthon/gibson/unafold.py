import os, subprocess, shlex, sys, csv



def self_prime(name, seq, t):
	cwd = os.getcwd()
	wd = '/home/bill/www/unafold/'
	os.chdir(wd)
	w = open(wd + name,'w')
	w.write(seq)
	w.close()
	devnull = open(os.devnull, 'w')
	cline = 'hybrid-ss-min -n DNA -t ' + t + ' -T ' + t + ' --mfold=5,-1,100 ' + wd + name
	p = subprocess.check_call(shlex.split(cline), stdout=devnull, stderr=devnull)
	cline = 'boxplot_ng ' + wd + name + '.plot'
	p = subprocess.check_call(shlex.split(cline), stdout=devnull, stderr=devnull)
	ss = csv.DictReader(open(wd + name + '.plot','r'), delimiter='\t')
	
	errors = []
	for r in ss:
		if r['j'] == str(len(seq)):
			errors.append((r['length'],float(r['energy'])/10))
	
	os.chdir(cwd)
	boxplot = 'woo'
	return errors, boxplot
	