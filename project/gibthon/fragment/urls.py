from django.conf.urls.defaults import patterns, include, url
from django.conf import settings

urlpatterns = patterns('fragment',
	(r'^$', 'views.fragments'),
	(r'^(\d+)/.*\.gb$', 'views.download'),
    (r'^(\d+)/.*$', 'views.fragment'),
    (r'^add$', 'views.add'),
	(r'^add/([A-Z]{2})$', 'views.add_submit'),
)

	