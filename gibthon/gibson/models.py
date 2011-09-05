# Gibson.models
# 
# Contains three significant classes:
#
# 1) Construct
# 	This is the main class of this app, containing all of the information 
#   about the construct.
#
# 2) ConstructFragment
#	This is an inbetween class linking each construct to its fragments. It
#	defines the start and end positions, and has methods to get the sequence
#	and list of features
#
# 3) Primer
#	This is the class for primers generated by the GCD
#
# Along with this there are a few further classes
#
# 1) PrimerHalf
#	Each primer is stored as two primerhalfs, one for the sticky end, one for the flappy end
#
# 2) Setings
# 	A class separate from the Construct class for storing settings for each 
#	construct. No huge benefit to it being separate from Construct, other than
#	separating fairly different aspects. Also makes forms a bit easier.
#
# 3) Warning
#	Applied to a primer when it is below optimal Tm, or has a mispriming
#
# And a couple of functions
#
# 1) hybrid_options
#	Generates the command line for using hybrid/hybrid-ss
#
# 2) add_fragment
#	Adds a new construct fragment to a construct. Should probably be a method of Construct

from django.db import models
from django import forms
from Bio.SeqUtils.MeltingTemp import Tm_staluc
from Bio.Seq import reverse_complement, Seq
from django.conf import settings


from Bio import SeqIO, Entrez
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from annoying.fields import AutoOneToOneField

from fragment.models import Gene

import os, subprocess, shlex, sys, csv
from subprocess import CalledProcessError

def hybrid_options(t,settings):
	return ' -n DNA -t %.2f -T %.2f -N %.2f -M %.2f --mfold=5,-1,100 ' %(t + settings.ss_safety, t + settings.ss_safety, settings.na_salt, settings.mg_salt)

def run_subprocess(cline, primer):
	try:
		p = subprocess.Popen(shlex.split(cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	except IOError:
		return False
	stdout, stderr = p.communicate()
	if p.returncode > 0:
		w = Warning.objects.create(
			primer = primer,
			type = 'sy',
			text = 'Primer calculation failed with error %d\n%s\n%s'%(p.returncode,stderr,stdout),
		)
		return False
	dir(p)
	return True

class Settings(models.Model):
	# automatically create a Settings object whne you make a new construct object
	construct = AutoOneToOneField('Construct', related_name='settings')
	# the below are various parameters explained in Gibson docs
	mg_salt = models.DecimalField(max_digits=3, decimal_places=2, default=0.0)
	na_salt = models.DecimalField(max_digits=3, decimal_places=2, default=1)
	ss_safety = models.PositiveSmallIntegerField(default=3)
	min_anneal_tm = models.PositiveSmallIntegerField(default=50)
	min_primer_tm = models.PositiveSmallIntegerField(default=60)
	min_overlap = models.PositiveSmallIntegerField(default=20)
	
	def __unicode__(self):
		return 'Settings for ' + self.construct.name

class PCRSettings(models.Model):
	construct = AutoOneToOneField('Construct', related_name='pcrsettings')
	repeats = models.PositiveSmallIntegerField(default=1)
	volume_each = models.DecimalField(max_digits=3, decimal_places=1, default=12.5)
	error_margin = models.PositiveSmallIntegerField(default=10)
	buffer_s = models.DecimalField(max_digits=3, decimal_places=1, default=10)
	buffer_d = models.DecimalField(max_digits=3, decimal_places=1, default=1)
	dntp_s = models.DecimalField(max_digits=3, decimal_places=1, default=10)
	dntp_d = models.DecimalField(max_digits=3, decimal_places=1, default=0.8)
	enzyme_s = models.DecimalField(max_digits=3, decimal_places=1, default=2.5)
	enzyme_d = models.DecimalField(max_digits=3, decimal_places=1, default=2.5)
	primer_d = models.DecimalField(max_digits=3, decimal_places=1, default=0.4)
	template_d = models.DecimalField(max_digits=4, decimal_places=1, default=100)
	
	def m(self):
		return self.repeats * (1+(self.error_margin/100))
	
	def buffer_v(self):
		return self.m() * self.volume_each * self.buffer_d/self.buffer_s
	
	def dntp_v(self):
		return self.m() * self.volume_each * self.dntp_d/self.dntp_s
	
	def enzyme_v(self):
		return self.m() * self.enzyme_d/self.enzyme_s
	
	def primer_v(self, primer_s):
		return self.m() * self.volume_each * self.primer_d/primer_s
	
	def template_v(self, template_s):
		return self.m() * self.template_d/template_s
	
	def water_v(self, primer_t_s, primer_b_s, template_s):
		return (self.m()*self.volume_each) - self.buffer_v() - self.dntp_v() - self.enzyme_v() - self.primer_v(primer_t_s) - self.primer_v(primer_b_s) - self.template_v(template_s)
	
	def total_v(self):
		return self.m() * self.volume_each
	

class Warning(models.Model):
	primer = models.ForeignKey('Primer', related_name='warning')
	# five warning types at the moment
	WARNING_TYPE = (
		('mp', 'MISPRIME'),
		('sp','SELF PRIME'),
		('ta','ANNEAL TM'),
		('tp','PRIMER TM'),
		('sy','SYSTEM ERROR'),
	)
	type = models.CharField(max_length=2, choices=WARNING_TYPE)
	text = models.CharField(max_length=150)

class Primer(models.Model):
	name = models.CharField(max_length=80)
	construct = models.ForeignKey('Construct', related_name='primer')
	flap = models.OneToOneField('PrimerHalf', related_name='flap')
	stick = models.OneToOneField('PrimerHalf', related_name='stick')
	boxplot = models.ImageField(upload_to='boxplots')
	concentration = models.DecimalField(default=5, max_digits=4, decimal_places=1)
	
	def vol(self):
		return self.construct.pcrsettings.primer_v(primer_s=self.concentration)
	
	class Meta:
		ordering = ['stick']
	
	def __unicode__(self):
		return self.name + ' (' + str(len(self.warning.all())) + ')'
		
	def csv(self):
		# returns an array of info to generate a .csv file
		return [self.name, self.length(), self.tm(), self.seq()]
	
	def length(self):
		return len(self.seq())
		
	def seq(self):
		return self.flap.seq() + self.stick.seq()
	
	def seq_pretty(self):
		# lowerUPPER case for printing the primer
		return str(self.flap.seq()).lower() + str(self.stick.seq()).upper()
		
	def tm(self):
		return round(Tm_staluc(str(self.seq())),2)
		
	def tm_len_anneal(self, target):
		# extends the length of the annealing portion of the primer until its target tm is reached
		while self.stick.tm() < target:
			if not self.stick.extend():
				w = Warning.objects.create(
					primer = self,
					type = 'tp',
					text = 'Anneal tm below target of ' + str(target) + '&deg;C ('+str(self.stick.tm())+'&deg;C)',
				)
				break
		self.stick.save()

	def tm_len_primer(self, target):
		# extends teh length of the flappy end of the primer until its target is reached
		while self.tm() < target:
			if not self.flap.extend():
				w = Warning.objects.create(
					primer = self,
					type = 'tp',
					text = 'Primer tm below target of ' + str(target) + '&deg;C ('+str(self.tm())+'&deg;C)',
				)
				break
		self.flap.save()
	
	def del_all(self):
		# delete the primer. requried because deleting self does not delete primerhalfs for some reason
		self.flap.delete()
		self.stick.delete()
		self.delete()
	
	def self_prime_check(self):
		self.warning.filter(type='sy').delete()
		# check for self prime events
		name = str(self.id)
		cwd = os.getcwd()
		# perform all work here. this should be in settings at some point
		wd = settings.UNAFOLD_WD
		os.chdir(wd)
		# dump the primer sequence into a file, because hyrbid* only supports file input
		w = open(wd + name,'w')
		w.write(str(self.seq()))
		w.close()
		devnull = open(os.devnull, 'w')
		# hyrbidise!
		cline = settings.HYBRID_SS_MIN_PATH + hybrid_options(self.tm(), self.construct.settings) + wd + name
		if not run_subprocess(cline, self):
			return False
		# generate a pretty boxplot
		cline = settings.BOXPLOT_NG_PATH + ' -t "Energy Dotplot for ' + name + ' " ' + wd + name + '.plot'
		if not run_subprocess(cline, self):
			return False
		# go through the boxplot info and check for mispriming
		try:
			csvfile = open(wd + name + '.plot', 'r')
		except IOError as e:
			w = Warning.objects.create(
					primer = self,
					type = 'sy',
					text = 'Error making boxplot',
				)
			return False
		ss = csv.DictReader(csvfile, delimiter='\t')
		warnings = []
		for r in ss:
			if int(r['j']) == len(self.seq()):
				warnings.append((r['length'],float(r['energy'])/10))
		# conver the boxplot to a png
		cline = 'convert ' + wd + name + '.ps ' + wd + name + '.png'
		if not run_subprocess(cline, self):
			return False		
		# move it to the media directory
		os.rename(wd+name+'.png',settings.MEDIA_ROOT+'unafold/'+name+'.png')
		# clean up after yourself
		for f in os.listdir('.'):
			if os.path.isfile(f) and f.startswith(name):
				os.remove(f)
		os.chdir(cwd)
		# delete old warnings and update info
		self.boxplot = name + '.png'
		self.warning.all().filter(type='sp').delete()
		self.save()
		for warning in warnings:
			w = Warning.objects.create(
				primer = self,
				type = 'sp',
				text = 'Potential self-priming of 3\' end! Length: ' + str(warning[0]) + ', dG: ' + str(warning[1]),
			)
			
			
	def corrprime(self):
		primer_name = str(self.name)
		fragment_name = str(self.construct.name) + '-' + str(self.stick.cfragment.fragment.name)
		cwd = os.getcwd()
		wd = settings.UNAFOLD_WD
		os.chdir(wd)
		w = open(wd + primer_name,'w')
		w.write(str(self.seq()))
		w.close()
		w = open(wd + fragment_name,'w')
		w.write(str(reverse_complement(self.stick.seq())))
		w.close()
		devnull = open(os.devnull, 'w')
		cline = settings.HYBRID_MIN_PATH + hybrid_options(self.tm(), self.construct.settings) + wd + fragment_name + ' ' + wd + primer_name
		p = subprocess.check_call(shlex.split(cline), stdout=devnull, stderr=devnull)
		
	
	def misprime_check(self):
		primer_name = str(self.id)
		fragment_name = str(self.construct.id) + '-' + str(self.stick.cfragment.id)
		cwd = os.getcwd()
		wd = settings.UNAFOLD_WD
		os.chdir(wd)
		w = open(wd + primer_name,'w')
		w.write(str(self.seq()))
		w.close()
		w = open(wd + fragment_name,'w')
		if(self.stick.top):
			w.write(str(self.stick.cfragment.sequence()))
		else:
			w.write(str(reverse_complement(Seq(self.stick.cfragment.sequence()))))
		w.close()
		devnull = open(os.devnull, 'w')
		cline = settings.HYBRID_MIN_PATH + hybrid_options(self.tm(), self.construct.settings) + wd + fragment_name + ' ' + wd + primer_name
		p = subprocess.check_call(shlex.split(cline), stdout=devnull, stderr=devnull)
		ss = csv.DictReader(open(wd + fragment_name + '-' + primer_name + '.plot'), delimiter='\t')
		warnings = []
		for r in ss:
			# length of annealing
			l = int(r['length'])
			# bp index from 3' end that is annealing
			j = len(self.seq()) - (int(r['j']) - len(self.stick.cfragment.sequence()))
			# bp index from 5' end of fragment
			i = (int(r['i']) + l - 1)
			if (j == 0 and i == len(self.stick.cfragment.sequence())) or (l == 1 or l == len(self.stick.seq())-1):
				# that's the priming we wanted! store the energy for comparison
				continue
			else:
				warnings.append((l,j,i,float(r['energy'])/10))
		self.warning.all().filter(type='mp').delete()
		for warning in warnings:
			w = Warning.objects.create(
				primer = self,
				type = 'mp',
				text= 'Potentital mis-priming ' + (str(warning[1]) + ' bp from ' if warning[1] > 0 else ' of ') + '3\' end of primer at bp ' + str(warning[2]) + ', length ' + str(warning[0]) + ', energy ' + str(warning[3]),
			)
		for f in os.listdir('.'):
			if os.path.isfile(f) and f.startswith(self.construct.name):
				os.remove(f)
		os.chdir(cwd)


class PrimerHalf(models.Model):
	cfragment = models.ForeignKey('ConstructFragment', related_name='ph')
	top = models.BooleanField()
	length = models.PositiveSmallIntegerField()
	
	def start(self):
		if self.top ^ self.isflap():
			return self.cfragment.end() - self.length
		else:
			return self.cfragment.start()
	
	def end(self):
		if self.top ^ self.isflap():
			return self.cfragment.end()
		else:
			return self.cfragment.start() + self.length
	
	def isflap(self):
		try: self.flap
		except: return False
		else: return True
		
	def extend(self):
		self.length += 1
		if self.start() < (self.cfragment.start_feature.start - self.cfragment.start_offset) or self.end() > (self.cfragment.end_feature.end + self.cfragment.end_offset):
			self.length -= 1
			return False
		else:
			return True
	
	def seq_surround(self):
		s = Seq(self.cfragment.fragment.sequence)
		start = max(self.start()-20,0)
		end = min(self.end()+20,len(self.cfragment.fragment.sequence))
		s = s[start:end]
		f = [self.cfragment.start()+self.cfragment.start_offset-1 <= i < self.cfragment.end()-self.cfragment.end_offset for i in range(start, end)]
		p = [self.start()-1 <= i < self.end() for i in range(start, end)]
		if self.top: 
			s = reverse_complement(s)
			f.reverse()
			p.reverse()
		bases = zip(s,f,p)
		s = ''
		for b in bases:
			s += '<td class="'+('feature ' if b[1] else '')+('primer' if b[2] else '') +'">'+b[0]+'</td>'
		return s
	
	def seq(self):
		s = Seq(self.cfragment.fragment.sequence)
		s = s[self.start():self.end()]
		if self.top: s = reverse_complement(s)
		return s
		
	def tm(self):
		return round(Tm_staluc(str(self.seq())),2)
		
	def __unicode__(self):
		if self.isflap():
			return self.flap.name + ' (flap): ' + str(self.seq()) + ' (' + str(self.tm()) + ')'
		else:
			return self.stick.name + ' (stick): ' + str(self.seq()) + ' (' + str(self.tm()) + ')'

SHAPE_CHOICES = (
	('c', 'Circular'),
	('l', 'Linear'),
)
class Construct(models.Model):
	owner = models.ForeignKey('auth.User', null=True)
	name = models.CharField(max_length=80)
	description = models.CharField(max_length=2000)
	genbank = models.OneToOneField('fragment.Gene', blank=True, null=True, related_name='construct_master')
	fragments = models.ManyToManyField('fragment.Gene', through='ConstructFragment', blank=True, related_name='construct_slave')

	shape = models.CharField(max_length=1, choices=SHAPE_CHOICES)
	created = models.DateTimeField(auto_now_add=True)
	modified = models.DateTimeField(auto_now=True)

	def __unicode__(self):
		return self.name
		
	def sequence(self):
		dna = ''
		for f in self.cf.all():
			dna += f.sequence()
		return dna
	
	def length(self):
		return len(self.sequence())
	
	def features(self):
		acc = 0
		for fr in self.cf.all():
			fe = fr.features()
			if fr.direction == 'r':
				fe.reverse()
			for f in fe:
				if fr.direction == 'r':
					t  = f.start
					f.start = fr.fragment.length() - f.end + 1
					f.end = fr.fragment.length() - t + 1
				f.start -= fr.start() - acc - 1
				f.end -= fr.start() - acc - 1
				yield f
			acc += fr.end() - fr.start() + 1
	
	def feature_count(self):
		return sum(1 for f in self.features())
	
	def features_pretty(self):
		acc = 0
		for fr in self.cf.all():
			fe = fr.features()
			if fr.direction == 'r':
				fe.reverse()
			for f in fe:
				if fr.direction == 'r':
					t  = f.start
					f.start = fr.fragment.length() - f.end + 1
					f.end = fr.fragment.length() - t + 1
				f.start -= fr.start() - acc - 1
				f.end -= fr.start() - acc - 1
				yield f.pretty() + (' [reverse]' if fr.direction == 'r' else '')
			acc += fr.end() - fr.start() + 1
	
	def gb(self):
		g = SeqRecord(
			Seq(self.sequence(),IUPAC.IUPACUnambiguousDNA()),
			id=self.name,
			name=self.name,
			description=self.description
		)
		g.features = [SeqFeature(FeatureLocation(ExactPosition(f.start-1),ExactPosition(f.end)), f.type, qualifiers=dict([[q.name,q.data] for q in f.qualifiers.all()])) for f in self.features()]
		return g.format('genbank')
		
	def add_fragment(self, fragment):
		o = len(self.fragments.all())
		cf = ConstructFragment.objects.create(construct=self, fragment=fragment, order = o, direction='f', start_feature=fragment.features.all()[0], end_feature=fragment.features.all()[0], start_offset=0, end_offset=0)
		if cf:
			return True
		else:
			return False
		
	def process(self, reset=True, new=True):
		if new:
			# delete all existing primers
			for p in self.primer.all():
				p.del_all()
		# used for returning progress
		n = self.cf.count()
		# reset offsets to zero
		if reset:
			for cf in self.cf.all():
				cf.start_offset = 0
				cf.end_offset = 0
				cf.save()
		if new:
			for i,cf in enumerate(self.cf.all()):
				cfu = self.cf.all()[(i+1)%n]
				pt = Primer.objects.create(
					name = self.name + '-' + cf.fragment.name + '-top',
					construct = self,
					stick = PrimerHalf.objects.create(
						cfragment = cf,
						top = True,
						length = self.settings.min_overlap
					),
					flap = PrimerHalf.objects.create(
						cfragment = cfu,
						top = True,
						length = self.settings.min_overlap
					)
				)
				cfd = self.cf.all()[(i-1)%n]
				pb = Primer.objects.create(
					name = self.name + '-' + cf.fragment.name + '-bottom',
					construct = self,
					stick = PrimerHalf.objects.create(
						cfragment = cf,
						top = False,
						length = self.settings.min_overlap
					),
					flap = PrimerHalf.objects.create(
						cfragment = cfd,
						top = False,
						length = self.settings.min_overlap
					)
				)
		else:
			for p in self.primer.all():
				p.stick.length = self.settings.min_overlap
				p.flap.length = self.settings.min_overlap
				p.save()
		for i,p in enumerate(self.primer.all()):
			if self.settings.min_anneal_tm > 0:
				p.tm_len_anneal(self.settings.min_anneal_tm)
			if self.settings.min_primer_tm > 0:
				p.tm_len_primer(self.settings.min_primer_tm)
			p.self_prime_check()
			yield ':%d'%(((2*i)+1)*(90.0/(4.0*n)))
			yield ' '*1024
			p.misprime_check()
			yield ':%d'%(((2*i)+2)*(90.0/(4.0*n)))
			yield ' '*1024
		yield ':100'		

class ConstructFragment(models.Model):
	construct = models.ForeignKey('Construct', related_name='cf')
	fragment = models.ForeignKey('fragment.Gene', related_name='cf')
	order = models.PositiveIntegerField()
	DIRECTION_CHOICES = (
		('f', 'Forward'),
		('r', 'Reverse'),
	)
	direction = models.CharField(max_length=1, choices=DIRECTION_CHOICES)
	start_feature = models.ForeignKey('fragment.Feature', related_name='start_feature')
	start_offset = models.IntegerField()
	end_feature = models.ForeignKey('fragment.Feature', related_name='end_feature')
	end_offset = models.IntegerField()
	concentration = models.DecimalField(default=100, max_digits=4, decimal_places=1)

	class Meta:
		ordering = ['order']
	
	def primer_top(self):
		for ph in self.ph.all():
			try: ph.stick
			except: continue
			else:
				if ph.top:
					return ph.stick
	
	def primer_bottom(self):
		for ph in self.ph.all():
			try: ph.stick
			except: continue
			else:
				if not ph.top:
					return ph.stick
	
	def features(self):
		feat = []
		for f in self.fragment.features.all():
			if self.direction == 'f':
				if f.start >= self.start() and f.end <=self.end():
						feat.append(f)
			else:
				if (self.fragment.length() - f.end + 1) >= self.start() and (self.fragment.length() - f.start + 1) <=self.end():
					feat.append(f)
		return feat
	
	def start(self):
		if self.direction == 'f':
			return self.start_feature.start - self.start_offset
		else:
			return self.fragment.length() - (self.start_feature.end - self.start_offset) + 1
	
	def end(self):
		if self.direction == 'f':
			return self.end_feature.end + self.end_offset
		else:
			return self.fragment.length() - (self.end_feature.start + self.end_offset) + 1
	
	def sequence(self):
		seq = self.fragment.sequence
		if self.direction == 'r':
			seq = str(reverse_complement(Seq(seq)))
		return seq[self.start()-1:self.end()]
	
	def tm(self):
		return ((self.primer_top().stick.tm() + self.primer_bottom().stick.tm())/2)-4
	
	def time(self):
		return (self.end()-self.start()+1)*45.0/1000
	
	def vol(self):
		return self.construct.pcrsettings.template_v(self.concentration)
	
	def water_v(self):
		return self.construct.pcrsettings.water_v(self.primer_top().concentration, self.primer_bottom().concentration, self.concentration)
	
	def __unicode__(self):
		return self.construct.name + ' : ' + self.fragment.name + ' (' + str(self.order) + ')'
	
