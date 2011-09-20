from django.conf.urls.defaults import patterns, include, url
from django.conf import settings

urlpatterns = patterns('molcal.views',
	(r'^$', 'molcalc'),
	(r'^get$', 'getall'),
)

