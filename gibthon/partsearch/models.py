from django.db import models

class PartType(models.Model):
	"""
	Every valid type of part
	- name: the name of the type. could this be the primary key also?		
	"""
	name = models.CharField(max_length=24, unique=True)

class Category(models.Model):
	"""
	Categories are handled as a directory structure
	-name: the full name of the Category (i.e. '//one/two/three')
	"""
	name = models.CharField(max_length=24)
#	parent = models.ForeignKey('self', related_name='children', null=True, blank=True)	
	
	def as_list(self):
		return self.name.strip('/').split('/')
	
	class Meta:
		verbose_name_plural = 'Categories'

class Part(models.Model):
	"""
	Represent a part, and it's relationships to other parts
	"""
	name = models.CharField("Name", max_length=16, unique=True)
	short_name = models.CharField("Short Name", max_length=12)
	short_desc = models.CharField("Short Description", max_length=512)
	type = models.ForeignKey(PartType)
	status = models.CharField(max_length = 24, blank=True)
	results = models.CharField(max_length = 24, blank=True)
	nickname = models.CharField(max_length = 24, blank=True)
	rating = models.DecimalField(blank=True, max_digits=5, decimal_places=3)
	url = models.URLField()
	date_entered = models.DateField()
	author = models.CharField(max_length = 256, blank=True)
	best_quality = models.CharField(max_length=32, blank=True)
	category = models.ForeignKey(Category)
	#TODO, write your own sequence object!
	sequence = models.CharField(max_length = 20000)
	
	deep_subparts = models.ManyToManyField('self', related_name="deep_superparts")
	specified_subparts = models.ManyToManyField('self', related_name="specified_subparts")
	twins = models.ManyToManyField('self', related_name="+")#can't go backwards
	
class Feature(models.Model):
	DIRECTIONS = (
		(1 , u'Forward'),
		(-1, u'Reverse'),
		(0 , u'Both'),
	)
	part = models.ForeignKey(Part)
	title = models.CharField(max_length=24)
	type = models.CharField(max_length=32)
	direction = models.IntegerField(choices = DIRECTIONS)
	#perhaps one day a location field to cope with genbank's fuzzyness?
	start = models.IntegerField()
	end = models.IntegerField()