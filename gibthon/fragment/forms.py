from django import forms

########################################## fields

class SequenceField(forms.CharField):
	def __init__(self, **kw):
		super(SequenceField, self).__init__(**kw)
		
	def clean(self, value):
		#clean the chars
		super(SequenceField, self).clean(value)
		#check id the sequence makes sense
		#assume IUPAC DNA
		from Bio.Seq import Seq
		from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
		letters = IUPACAmbiguousDNA().letters.lower()
		if not letters: 
			raise ValueError("Alphabet does not define letters.") 
		for letter in value.lower(): 
			if letter not in letters: 
				raise forms.ValidationError(u"Invalid sequence, IUPAC Ambiguous DNA only accepts '%s'" % letters)
		#sequence ok
		return Seq(value, IUPACAmbiguousDNA()) 


############################################# forms

class MetaForm(forms.Form):
	name = forms.CharField(	max_length=32, 
							widget=forms.widgets.TextInput(attrs={'cols': 32,})
							)
	desc = forms.CharField(	required = False,
							max_length=4096, 
							label="Description", 
							widget=forms.widgets.Textarea(
								attrs={'cols': 60, 'rows':4,}))
	
class SequenceForm(forms.Form):
	seq = SequenceField(	max_length=500000, 
							label="Sequence",
							widget=forms.widgets.Textarea(
								attrs={'cols': 60,'rows':11,})
							)

class FragmentForm(MetaForm, SequenceForm):
	pass

