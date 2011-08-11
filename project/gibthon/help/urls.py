from django.conf.urls.defaults import patterns, include, url
from django.conf import settings

urlpatterns = patterns('help.views',
	(r'^$', 'help'),
	(r'^gibson/$', 'gibson'),
	(r'^primer/$', 'primer'),
	(r'^gibthon/$', 'gibthon'),
)