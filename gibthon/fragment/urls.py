from django.conf.urls.defaults import patterns, include, url
from django.conf import settings

from importfragment import entrezpatterns

urlpatterns = patterns('fragment',
	(r'^$', 'views.fragments'),
	(r'^(\d+)/.*\.gb$', 'views.download'),
	(r'^download/$', 'views.download_multi'),
	(r'^get/(\d+)/$', 'api.get'),
	(r'^(\d+)/.*$', 'views.fragment'),
	(r'^(\d+)/$', 'views.fragment'),
	(r'^add$', 'views.add'),
	(r'^delete/$', 'views.delete'),
	(r'^import/$', 'importfragment.fragment_import'),
	(r'^import/entrez/', include(entrezpatterns)),
	
)

	
