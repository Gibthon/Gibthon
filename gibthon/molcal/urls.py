from django.conf.urls.defaults import patterns, include, url
from django.conf import settings

urlpatterns = patterns('molcal.views',
	(r'^$', 'molcal'),
	(r'^get$', 'getall'),
)

