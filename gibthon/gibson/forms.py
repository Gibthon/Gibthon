# gibson.forms
#
# contains forms for use with the gibthon app, along with classes for radio
# buttons that are compatible with jQuery buttonsets. These should probably be
# moved somewhere a bit more global at some point

from django import forms
from models import *

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

# very basic form for changing settings
class SettingsForm(forms.ModelForm):
	class Meta:
		model = Settings
		exclude = ['construct']

# for changing the detail of the form
class ConstructForm(forms.ModelForm):
	description = forms.CharField(widget=forms.Textarea)
	shape = forms.ChoiceField(
		widget=forms.RadioSelect(renderer = BetterRadioFieldRenderer),
		choices=SHAPE_CHOICES,
		initial='c'
	)
	class Meta:
		model = Construct
		exclude = ['genbank', 'fragments', 'settings', 'owner']
		
# for generating the content of the accordion used to manipulate ConstructFragments
class FeatureListForm(forms.Form):
	DIRECTION_CHOICES = (
		('f', 'Forward'),
		('r', 'Reverse'),
	)
	start_feature = forms.ModelChoiceField('fragment.Feature', None, label='')
	finish_feature = forms.ModelChoiceField('fragment.Feature', None, label='')
	direction = forms.ChoiceField(widget=forms.RadioSelect(renderer = BetterRadioFieldRenderer), choices=DIRECTION_CHOICES)
	def __init__(self, _fragment, _construct, *args, **kwargs):
		sf = self.base_fields['start_feature']
		ff = self.base_fields['finish_feature']
		cf = _fragment.cf.get(fragment=_fragment.pk, construct=_construct.pk)
		sf.queryset = _fragment.features.all()
		ff.queryset = _fragment.features.all()
		sf.widget.choices = sf.choices
		ff.widget.choices = ff.choices
		sf.initial = cf.start_feature
		ff.initial = cf.end_feature
		print cf.direction == 'r'
		self.base_fields['direction'].initial = cf.direction
		super(FeatureListForm, self).__init__(*args, **kwargs)
