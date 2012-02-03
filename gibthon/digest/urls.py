from django.conf.urls.defaults import patterns, include, url
from django.conf import settings

urlpatterns = patterns('digest.views',
	( r'^$', 'digest' ),
)
