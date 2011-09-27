from Bio import SeqIO
from Bio.Alphabet import IUPAC

import os, subprocess, shlex, sys, csv, random, string
from subprocess import CalledProcessError


class UnaFolder():
	def __init__(self, una_root = '/usr/local/bin/', wd = '/tmp/', t=60, safety=3, na_salt=1, mg_salt=0):
		self.una_root = una_root
		self.wd = wd
		self.name = ''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for x in range(20))
		self.t = t
		self.safety = safety
		self.na_salt = na_salt
		self.mg_salt = mg_salt
		
	def hybrid_options(self):
		return ' -n DNA -t %.2f -T %.2f -N %.2f -M %.2f --mfold=5,-1,100 ' %(self.t + self.safety, self.t + self.safety, self.na_salt, self.mg_salt)
	
	def process(self, cline):
		try:
			p = subprocess.Popen(shlex.split(cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		except IOError as e:
			print e
			return False
		stdout, stderr = p.communicate()
		dir(p)
		if p.returncode > 0:
			print("oh return code")
			return False
		return True
	
	def self_prime(self, sequence):
		cwd = os.getcwd()
		os.chdir(self.wd)
		w = open(self.wd + self.name, 'w')
		w.write(str(sequence))
		w.close()
		devnull = open(os.devnull,'w')
		cline = self.una_root + 'hybrid-ss-min' + self.hybrid_options() + self.wd + self.name
		if not self.process(cline):
			print('Could not hybridise')
			return False
		cline = self.una_root + 'boxplot_ng -t "Energy dotplot "' + self.wd + self.name + '.plot'
		if not self.process(cline):
			print('Could not convert boxplot')
			return False
		try:
			csvfile = open(self.wd + self.name + '.plot', 'r')
		except IOError as e:
			print e
			return False
		ss = csv.DictReader(csvfile, delimiter='\t')
		warnings = []
		for r in ss:
			if int(r['j']) == len(self.seq()):
				warnings.append(r['length'], float(r['energy']/10))
		cline = 'convert ' + self.wd + self.name + '.ps ' + self.wd + self.name + '.png'
		if not self.process(cline):
			return False
		for f in os.listdir(self.wd):
			if os.path.isfile(f) and f.startswith(self.name):
				os.remove(f)
		os.chdir(cwd)