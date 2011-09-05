from django.db import models
from django import forms
from django.core.exceptions import ObjectDoesNotExist


from Bio import SeqIO, Entrez
from Bio import SeqFeature
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
	sequence = models.TextField(max_length=500000)

	ORIGIN_CHOICES = (
		('NT', 'Nucleotide Database'),
		('BB', 'BioBrick'),
		('UL', 'Upload'),
		('GC', 'Gibthon Construct'),
		('MN', 'Manual'),
	)
	origin = models.CharField(max_length=2, choices=ORIGIN_CHOICES)

	def __unicode__(self):
		return self.name
	

	def add(_record, _origin, _user, errors = False):
		t = 0
		g = Gene(owner=_user, name=_record.name,description=_record.description,sequence=_record.seq, origin=_origin)
		g.save()
		for key,value in _record.annotations.items():
			try:
				Annotation.add(g, key, value)
			except Exception as ex:
				if ex.message.lower().startswith('data truncated'):
					t += 1
				else:
					print "Error: %s" % ex.message 
		for feature in _record.features:
			Feature.add(feature,g,_origin)
		
		if(errors):
			e = {'truncated': t,}
			return (g, e)
		
		return g
	add = staticmethod(add)
	
	def remove(_owner, _id):
		try:
			g = Gene.objects.get(owner=_owner, pk=_id)
			Annotation.remove(g)
			Feature.remove(g)
			g.delete()
			return True
		except ObjectDoesNotExist:
			return False
	remove = staticmethod(remove)
	
	def gb(self):
		g = SeqRecord(Seq(self.sequence,IUPAC.IUPACUnambiguousDNA()),id=self.name, name=self.name, description=self.description)
		g.features = [SeqFeature(FeatureLocation(ExactPosition(f.start-1),ExactPosition(f.end)), f.type, qualifiers=dict([[q.name,q.data] for q in f.qualifiers.all()])) for f in self.features.all()]
		return g.format('genbank')
	
	def to_seq_record(self):
		"""Convert the Gene to a SeqRecord"""
		#build a list of features
		feats = [_f.to_seq_feature() for _f in self.features.all()]
		#build a dictionary of annotations & refs
		annot = {}
		for a in self.annotations.all():
			a.to_ann(annot)
		annot['references'] = [r.to_ref() for r in self.references.all()]
		
		return SeqRecord(	seq = Seq(self.sequence, IUPAC.IUPACAmbiguousDNA()),
								name = self.name,
								description = self.description,
								features = feats,
								annotations = annot)		
		
	def pretty(self):
		seg = 8
		segline = 10
		s = [self.sequence[j:j+(segline*seg)] for j in range(0, len(self.sequence), segline*seg)]
		i = range(1, len(self.sequence), segline*seg)
		s = [SeqLine(ii,ss) for ss,ii in zip(s,i)]
		for f in self.features.all():
			i = 0
			while s[i].number < f.start - (seg*segline):
				i += 1
			j=i
			while s[j].number < f.end - (seg*segline):
				j += 1
			q = f.qualifiers.all()[0].data if f.qualifiers.all() else ''
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
	
	def to_ref(self):
		r = SeqFeature.Reference()
		r.title = self.title
		r.authors = self.authors
		r.journal = self.journal
		r.medline_id = self.medline_id
		r.pubmed_id = self.pubmed_id
		return r
	
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
				a = Annotation(gene= _gene, key=_key, value = str(v))
				a.save()
		else:
			a = Annotation(gene = _gene, key =str(_key), value = str(_value))
			a.save()
	add = staticmethod(add)
	
	def to_ann(self, d):
		"""add the annotation to d"""
		d[self.key] = self.value
	
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
	start = models.PositiveIntegerField() #This is stored as 0-offset
	end = models.PositiveIntegerField() #This is stored as 0-offset
	DIRECTION_CHOICES = (
		('f', 'Forward'),
		('r', 'Reverse'),
	)
	direction = models.CharField(max_length=1, choices=DIRECTION_CHOICES)
	gene = models.ForeignKey('Gene', related_name="features")
	
	class Meta:
		ordering = ['start']
	
	def add(feature, g, origin):
		#if (origin == "BB"): #but parts are 1-offset based!
		#	feature.location.start.position += 1
		if (feature.location.end.position == 0):
			return
		f = Feature(type=feature.type, start=feature.location.start.position, end=feature.location.end.position,
			direction='f' if feature.location.start < feature.location.end else 'r', gene = g)
		f.save()
		for _name,_data in feature.qualifiers.iteritems():
			data = ''
			if hasattr(_data, '__iter__'):
				data= ', '.join(_data)
			else:
				data=str(_data)
			q = Qualifier(name=_name,data=data,feature = f)
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
	
	def to_seq_feature(self):
		quals = {}
		for q in self.qualifiers.all():
			quals[q.name] = q.data
		s = None
		if self.direction == 'f':
			s = 1
		elif self.direction == 'r':
			s = -1
		return SeqFeature.SeqFeature(	location = SeqFeature.FeatureLocation(self.start, self.end),
												type = self.type,
												strand = s,
												qualifiers = quals) 
	
	def __unicode__(self):
		pos = ' (' + str(self.start) + '..' + str(self.end) + ')'
		if self.qualifiers.all():
			return self.qualifiers.all()[0].data + pos
		else:
			return self.type + pos
	
	def pretty(self):
		pos = ' (' + str(self.start) + '..' + str(self.end) + ')'
		if self.qualifiers.all():
			return self.qualifiers.all()[0].data + pos
		else:
			return self.type + pos
		
	def first(self):
		return self.qualifiers.all()[0]
	
class Qualifier(models.Model):
	name = models.CharField(max_length=30)
	data = models.CharField(max_length=512)
	feature = models.ForeignKey('Feature', related_name="qualifiers")
	
	def __unicode__(self):
		return self.name

	
