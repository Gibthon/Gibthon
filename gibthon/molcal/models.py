from django.db import models

# Create your models here.

class Mol(models.Model):
	name = models.CharField(max_length=80)
	formula = models.CharField(max_length=80)
	fullname = models.CharField(max_length=300)
	weight = models.DecimalField(max_digits=8, decimal_places=3)
	
	def __unicode__(self):
		return self.name