from django.conf.urls.defaults import patterns, include, url
from django.conf import settings

urlpatterns = patterns('pricalc.views',
	(r'^$', 'pricalc'),
	(r'^go$', 'go'),
)

