from django.db import models
from Bio.Seq import reverse_complement, Seq
from Bio.SeqUtils.MeltingTemp import Tm_staluc
from gibthon.unafold import UnaFolder
from django.conf import settings

import os

class Oligo():
	def __init__(self, data):
		self.left=data.get('gene1')
		self.right=data.get('gene2')
		self.sl=int(data.get('sl'))-1
		self.sr=int(data.get('sr'))-1
		self.el=int(data.get('el'))
		self.er=int(data.get('er'))
	
	def leftHalf(self):
		return self.left[::-1][self.sl:self.el][::-1]
	
	def rightHalf(self):
		return self.right[self.sr:self.er]
	
	def topPrimer(self):
		return str(reverse_complement(Seq(self.bottomPrimer())))
	
	def bottomPrimer(self):
		return self.leftHalf() + self.rightHalf()
		
	def topTm(self):
		return Tm_staluc(self.leftHalf())
		
	def bottomTm(self):
		return Tm_staluc(self.rightHalf())
	
	def fullTm(self):
		return Tm_staluc(self.bottomPrimer())
	
	def selfPrime(self):
		u = UnaFolder(t = self.fullTm(), safety=3, mg_salt=0, na_salt=1)
		ret, image = u.self_prime(str(self.topPrimer()))
		print settings.MEDIA_ROOT
		os.rename(image, settings.MEDIA_ROOT+image[1::])
		wst = ''
		for warning in u.warnings:
			wst += 'Potential self-priming of 3\' end! Length: ' + str(warning[0]) + ', dG: ' + str(warning[1]) + '<br />'
		return settings.MEDIA_URL + image, wst
	
		