# gibson.forms
#
# contains forms for use with the gibthon app, along with classes for radio
# buttons that are compatible with jQuery buttonsets. These should probably be
# moved somewhere a bit more global at some point

from django import forms
from models import *
from gibthon import formfields

# very basic form for changing settings
class SettingsForm(forms.ModelForm):
	class Meta:
		model = Settings
		exclude = ['construct']

# for changing the detail of the form
class ConstructForm(forms.ModelForm):
	description = forms.CharField(widget=forms.Textarea)
	shape = forms.ChoiceField(
		widget=forms.RadioSelect(renderer = formfields.BetterRadioFieldRenderer),
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
	direction = forms.ChoiceField(widget=forms.RadioSelect(renderer = formfields.BetterRadioFieldRenderer), choices=DIRECTION_CHOICES)
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
		self.base_fields['direction'].initial = cf.direction
		super(FeatureListForm, self).__init__(*args, **kwargs)
