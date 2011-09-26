from django.conf.urls.defaults import patterns, include, url
from django.conf import settings

importpatterns = patterns('gibthon.fragment.importfragment',
	(r'^entrez/$', 'entrez'),
	(r'^entrez/search/$', 'entrez_search'),
	(r'^entrez/summary/$', 'entrez_summary'),
	(r'^entrez/import/$', 'entrez_import'),
	(r'^upload/$', 'upload_form'),
	(r'^upload/go/$', 'handle_upload'),
	(r'^part/$', 'part_form'),
	(r'^part/go/$', 'part_import'),
	(r'^manual/$', 'manual_form'),
	(r'^manual/add/$', 'manual_add'),
)


urlpatterns = patterns('fragment',
	(r'^$', 'views.fragments'),
	(r'^(\d+)/.*\.gb$', 'views.download'),
	(r'^download/$', 'views.download_multi'),
	(r'^get/(\d+)/$', 'api.get'),
	(r'^list/$', 'api.get_fragments'),
	(r'^(\d+)/getMeta/$', 'views.fragment_meta'),
	(r'^(\d+)/getSeq/$', 'views.fragment_seq'),
	(r'^(\d+)/saveMeta/$', 'api.save_meta'),
	(r'^(\d+)/$', 'views.fragment'),
	(r'^delete/$', 'views.delete'),
	(r'^import/$', 'importfragment.fragment_import'),
	(r'^import/', include(importpatterns)),	
)

	
