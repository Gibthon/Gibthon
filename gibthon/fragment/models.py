from django.db import models
from django import forms
from django.core.exceptions import ObjectDoesNotExist


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
	sequence = models.TextField(max_length=50000)
	ORIGIN_CHOICES = (
		('NT', 'Nucleotide Database'),
		('BB', 'BioBrick'),
		('UL', 'Upload'),
		('GC', 'Gibthon Construct'),
	)
	origin = models.CharField(max_length=2, choices=ORIGIN_CHOICES)

	def __unicode__(self):
		return self.name
	
	def add(_record, _origin, _user):
		g = Gene(owner=_user, name=_record.name,description=_record.description,sequence=_record.seq, origin=_origin)
		g.save()
		for key,value in _record.annotations.items():
			Annotation.add(g, key, value)
		for feature in _record.features:
			f = Feature.add(feature,g,_origin)
		return g
	add = staticmethod(add)
	
	def remove(_owner, _id):
		try:
			g = Gene.objects.get(owner=_owner, pk=_id)
			Annotation.remove(g)
			Feature.remove(g)
			g.delete()
		except ObjectDoesNotExist:
			pass
	remove = staticmethod(remove)
	
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

class Reference(models.Model):
	"""Store a reference"""
	gene = models.ForeignKey('Gene', related_name='references')
	title = models.CharField(max_length=1024)
	authors = models.CharField(max_length=1024)
	journal = models.CharField(max_length=512)
	medline_id = models.CharField(max_length=24, blank=True)
	pubmed_id = models.CharField(max_length=24, blank=True)
	
	def add(_gene, _refs):
		for r in _refs:
			print "Adding Ref '%s': journal '%s'" % (r.title, r.journal)
			ref = Reference(	gene = _gene, 
									title = str(r.title),
									journal = str(r.journal) )
			if isinstance(r.authors, list):
				ref.authors = [str(author) for author in r.authors]
			else:
				ref.authors = r.authors
			
			if hasattr(r, 'medline_id'):
				ref.medline_id = r.medline_id
			
			if hasattr(r, 'pubmed_id'):
				ref.pubmed_id = r.pubmed_id
			
			ref.save()
	add = staticmethod(add)
	
	def remove(_gene):
		"""remove all refs concerning _gene"""
		try:
			Reference.objects.filter(gene=_gene).delete()
		except ObjectDoesNotExist:
			pass
	remove = staticmethod(remove)

class Annotation(models.Model):
	gene = models.ForeignKey('Gene', related_name='annotations')
	key = models.CharField(max_length=30)
	value = models.CharField(max_length=5120, blank=True)

	def add(_gene, _key, _value):
		"""add an annotation. Does not accept references"""
		if _key.lower() == 'references':
			Reference.add(_gene, _value)
			return
		if hasattr(_value, '__iter__'):
			for v in _value:
				print "Adding multi annotation '%s': '%s'" % (_key, v)
				a = Annotation(gene= _gene, key=_key, value = str(v))
				a.save()
		else:
			print "Adding annotation '%s': '%s'" % (_key, _value)
			a = Annotation(gene = _gene, key =str(_key), value = str(_value))
			a.save()
	add = staticmethod(add)
	
	def remove(_gene):
		"""remove all annotations relating to gene"""
		Reference.remove(_gene)
		try:
			obj = Annotation.objects.filter(gene=_gene)
			obj.delete()
		except ObjectDoesNotExist:
			pass
	remove = staticmethod(remove)
		
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
	
	def remove(_gene):
		"""remove all features of _gene"""
		try:
			feats = Feature.objects.filter(gene=_gene)
			for f in feats:
				try:
					Qualifier.objects.filter(feature=f).delete()
				except ObjectDoesNotExist:
					pass
			feats.delete()
		except ObjectDoesNotExist:
			pass
	remove = staticmethod(remove)
	
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
	data = models.CharField(max_length=512)
	feature = models.ForeignKey('Feature', related_name="qualifier")
	
	def __unicode__(self):
		return self.name


class BBForm(forms.Form):
	id = forms.CharField(max_length=30)
	
class NTForm(forms.Form):
	id = forms.CharField(max_length=30)
	
class ULForm(forms.Form):
	file = forms.FileField()
	
