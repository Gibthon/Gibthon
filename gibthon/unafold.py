from Bio import SeqIO
from Bio.Alphabet import IUPAC

import os, subprocess, shlex, sys, csv, random, string
from subprocess import CalledProcessError

import settings

class UnaFolder():
	def __init__(self, wd = '/tmp/', t=60, safety=3, na_salt=1, mg_salt=0):
		self.wd = settings.UNAFOLD_WD
		self.name = ''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for x in range(20))
		self.t = t
		self.safety = safety
		self.na_salt = na_salt
		self.mg_salt = mg_salt
		self.warnings = []
		
	def hybrid_options(self):
		return ' -n DNA -t %.2f -T %.2f -N %.2f -M %.2f --mfold=5,-1,100 ' %(self.t + self.safety, self.t + self.safety, self.na_salt, self.mg_salt)
	
	def process(self, cline):
		try:
			print 'subprocess.Popen(shlex.split(%s), stdout=subprocess.PIPE, stderr=subprocess.PIPE)' % cline
			p = subprocess.Popen(shlex.split(cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		except IOError as e:
			print e
			return False
		stdout, stderr = p.communicate()
		if p.returncode > 0:
			return False
		return True
	
	def self_prime(self, sequence):
		cwd = os.getcwd()
		os.chdir(self.wd)
		w = open(self.wd + self.name, 'w')
		w.write(str(sequence))
		w.close()
		cline = settings.HYBRID_SS_MIN_PATH + self.hybrid_options() + self.wd + self.name
		if not self.process(cline):
			print('Could not hybridise')
			return (False, None)
		cline = '%s -t "Energy Dotplot " %s%s.plot' % (settings.BOXPLOT_NG_PATH, self.wd, self.name)
		if not self.process(cline):
			print('Could not convert boxplot')
			return (False, None)
		try:
			csvfile = open(self.wd + self.name + '.plot', 'r')
		except IOError as e:
			print e
			return (False, None)
		ss = csv.DictReader(csvfile, delimiter='\t')
		warnings = []
		for r in ss:
			if int(r['j']) == len(sequence):
				self.warnings.append((r['length'], float(r['energy'])/10))
		cline = 'convert ' + self.wd + self.name + '.ps ' + self.wd + self.name + '.png'
		if not self.process(cline):
			return (False, None)
		for f in os.listdir(self.wd):
			if os.path.isfile(f) and f.startswith(self.name) and not f.endswith('.png'):
				os.remove(f)
		os.chdir(cwd)
		return (True, self.wd + self.name + '.png')
		
	def mis_prime(self, target, primer):
		primerloc = self.name + '-primer'
		targetloc = self.name + '-target'
		cwd = os.getcwd()
		os.chdir(self.wd)
		w = open(self.wd + primerloc, 'w')
		w.write(str(primer))
		w.close
		w = open(self.wd + targetloc, 'w')
		w.write(str(target))
		w.close()
		cline = settings.HYBRID_MIN_PATH + self.hybrid_options() + self.wd + targetloc + ' ' + self.wd + primerloc
		if not self.process(cline):
			print('Could not hybridise')
			return False
		try:
			csvfile = open(self.wd + targetloc + '-' + primerloc + '.plot')
		except IOError as e:
			print e
			return False
		ss = csv.DictReader(csvfile, delimiter='\t')
		for r in ss:
			l = int(r['length'])
			j = len(primer) - (int(r['j']) - len(target))
			i = (int(r['i']) + l - 1)
			if (j == 0 and i == len(target)) or (l == 1 or l == len(primer)-1):
				continue
			else:
				self.warnings.append((l,j,i,float(r['energy'])/10))
		
		os.chdir(cwd)
		return True
