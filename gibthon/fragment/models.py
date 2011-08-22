from django.db import models
from django import forms

from Bio import SeqIO, Entrez
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import urllib

Entrez.email = 'entrez@gibthon.org'
Entrez.tool = 'gibthon/biopython'

class SeqLine:
	def __init__(self,_number,_seq):
		self.number = _number
		self.seq = Seq(_seq, IUPAC.IUPACUnambiguousDNA())
		self.rseq = self.seq.complement()
		self.features = []
	
class LineFeature:
	def __init__(self,_start,_finish,_type,_label=''):
		ll = 8 * 10
		if _start == 'start':
			blank = ''
			if _finish == 'end':
				star = '*'*(ll)
			else:
				star = '*'*(_finish%(ll))
		else:
			blank = '&nbsp;'*(_start%(ll)-1)
			if _finish == 'end':
				star = '*'*(ll-_start%(ll)+1)
			else:
				star = '*'*(_finish-_start+1)
		self.block = blank + star
		self.type = blank + _type
		if _label:
			self.type += '(' + _label + ')'

class Gene(models.Model):
	owner = models.ForeignKey('auth.User')
	name = models.CharField(max_length=100)
	description = models.CharField(max_length=500)
	sequence = models.TextField()
	ORIGIN_CHOICES = (
		('NT', 'Nucleotide Database'),
		('BB', 'BioBrick'),
		('UL', 'Upload'),
		('GC', 'Gibthon Construct'),
	)
	origin = models.CharField(max_length=2, choices=ORIGIN_CHOICES)

	def __unicode__(self):
		return self.name
	
	def add(_data, _origin, _user):
		if _origin == "BB":
			gbf = urllib.urlretrieve('http://www.cambridgeigem.org/gbdownload/'+_data+'.gb')[0]
		elif _origin == "NT":
			handle = Entrez.esearch(db='nucleotide', term=_data)
			record = Entrez.read(handle)
			gbid = record['IdList'][0]
			handle = Entrez.efetch(db="nucleotide", id=gbid, rettype="gb")
			gbf = handle.fp
		elif _origin == "UL":
			gbf = _data
		record = SeqIO.read(gbf,'genbank')
		g = Gene(owner=_user, name=record.name,description=record.description,sequence=record.seq, origin=_origin)
		if len(record.seq > 100000):
			return False
		g.save()
		for feature in record.features:
			f = Feature.add(feature,g,_origin)
		return True
	add = staticmethod(add)
	
	def gb(self):
		g = SeqRecord(Seq(self.sequence,IUPAC.IUPACUnambiguousDNA()),id=self.name, name=self.name, description=self.description)
		g.features = [SeqFeature(FeatureLocation(ExactPosition(f.start-1),ExactPosition(f.end)), f.type, qualifiers=dict([[q.name,q.data] for q in f.qualifier.all()])) for f in self.feature.all()]
		return g.format('genbank')
		
	def pretty(self):
		seg = 8
		segline = 10
		s = [self.sequence[j:j+(segline*seg)] for j in range(0, len(self.sequence), segline*seg)]
		i = range(1, len(self.sequence), segline*seg)
		s = [SeqLine(ii,ss) for ss,ii in zip(s,i)]
		for f in self.feature.all():
			i = 0
			while s[i].number < f.start - (seg*segline):
				i += 1
			j=i
			while s[j].number < f.end - (seg*segline):
				j += 1
			q = f.qualifier.all()[0].data if f.qualifier.all() else ''
			if i==j:
				s[i].features.append(LineFeature((f.start if f.start > 0 else 'start'),f.end,f.type,q))
			else:
				s[i].features.append(LineFeature(f.start,'end',f.type,q))
				i += 1
				while i<j:
					s[i].features.append(LineFeature('start','end',f.type))
					i += 1
				s[i].features.append(LineFeature('start',f.end,f.type,q))
		return s
	
	def length(self):
		return len(self.sequence)
		
	
class Feature(models.Model):
	type = models.CharField(max_length=30)
	start = models.PositiveIntegerField()
	end = models.PositiveIntegerField()
	DIRECTION_CHOICES = (
		('f', 'Forward'),
		('r', 'Reverse'),
	)
	direction = models.CharField(max_length=1, choices=DIRECTION_CHOICES)
	gene = models.ForeignKey('Gene', related_name="feature")
	
	def add(feature, g, origin):
		if (origin == "BB"):
			feature.location.start.position += 1
		if (feature.location.end.position == 0):
			return
		f = Feature(type=feature.type, start=feature.location.start.position, end=feature.location.end.position,
			direction='f' if feature.location.start < feature.location.end else 'r', gene = g)
		f.save()
		for _name,_data in feature.qualifiers.iteritems():
			q = Qualifier(name=_name,data=_data[0],feature = f)
			q.save()
		return f
	add = staticmethod(add)
	
	def __unicode__(self):
		pos = ' (' + str(self.start) + '..' + str(self.end) + ')'
		if self.qualifier.all():
			return self.qualifier.all()[0].data + pos
		else:
			return self.type + pos
	
	def pretty(self):
		pos = ' (' + str(self.start) + '..' + str(self.end) + ')'
		if self.qualifier.all():
			return self.qualifier.all()[0].data + pos
		else:
			return self.type + pos
		
	def first(self):
		return self.qualifier.all()[0]
	
class Qualifier(models.Model):
	name = models.CharField(max_length=30)
	data = models.CharField(max_length=150)
	feature = models.ForeignKey('Feature', related_name="qualifier")
	
	def __unicode__(self):
		return self.name

class BBForm(forms.Form):
	id = forms.CharField(max_length=30)
	
class NTForm(forms.Form):
	id = forms.CharField(max_length=30)
	
class ULForm(forms.Form):
	file = forms.FileField()
	
