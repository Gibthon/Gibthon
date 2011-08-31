# gibthon.formfields - custom form fields

from django import forms
from django.utils.safestring import mark_safe
from django.utils.encoding import force_unicode
from django.utils.html import conditional_escape

import random # for generating unique ids for radio buttons

# classes for the Radio buttons compatible with jQuery buttonsets

class BetterRadioFieldRenderer(forms.widgets.RadioFieldRenderer):

	def __init__(self, *args, **kwargs):
		super(BetterRadioFieldRenderer, self).__init__(*args, **kwargs)
	
	# on pages with multiple forms, required unique ID for each buttonset.
	# appending a random four diget number is a fairly shoddy but effective way
	# of doing this
	def __iter__(self):
		key = str(random.randint(1000,9999))
		for i, choice in enumerate(self.choices):
			yield BetterRadioInput(key, self.name, self.value, self.attrs.copy(), choice, i)

	def __getitem__(self, idx):
		choice = self.choices[idx] # Let the IndexError propogate
		return BetterRadioInput(self.name, self.value, self.attrs.copy(), choice, idx)
	
	def render(self):
		return mark_safe(u'\n'.join([u'%s' % force_unicode(w) for w in self]))

class BetterRadioInput(forms.widgets.RadioInput):
	def __init__(self, key, *args, **kwargs):
		super(BetterRadioInput, self).__init__(*args, **kwargs)
		self.name += key
		self.attrs['id'] += '_' + key
	def __unicode__(self):
		if 'id' in self.attrs:
			label_for = ' for="%s_%s"' % (self.attrs['id'], self.index)
		else:
			label_for = ''
		choice_label = conditional_escape(force_unicode(self.choice_label))
		return mark_safe(u'%s<label%s>%s</label>' % (self.tag(), label_for, choice_label))
