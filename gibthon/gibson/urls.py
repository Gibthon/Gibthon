from django.conf.urls.defaults import patterns, include, url
from django.conf import settings

primerpatterns = patterns('gibson.views',
	(r'^$', 'primers'),
	(r'^download$', 'primer_download'),
	(r'^pdf$', 'pdf'),
	(r'^save$', 'primer_save'),
	(r'^reset$', 'primer_reset'),
	(r'^csv$', 'csv_primers'),
	(r'^(?P<pid>\d+)/[\w\-_]+/$', 'primer'),
	(r'^(?P<pid>\d+)$', 'load_primer'),
	(r'^(?P<pid>\d+)/[\w\-_]+/offset$', 'primer_offset'),
	(r'^\d+/[\w\-_]+/reset$', 'primer_reset'),
)

constructpatterns = patterns('gibson.views',
	(r'^$', 'construct'),
	(r'^summary/$', 'summary'),
	(r'^settings/$', 'construct_settings'),
	(r'^delete$', 'construct_delete'),
	(r'^view/(?P<fid>\d+)$', 'fragment_viewer'),
	(r'^add$', 'fragment_browse'),
	(r'^add/(?P<fid>\d+)$', 'fragment_add'),
	(r'^delete/(?P<cfid>\d+)$', 'fragment_delete'),
	(r'^save/$', 'save'),
	(r'^process/$', 'process'),
	(r'^primers/', include(primerpatterns)),
)

designerpatterns = patterns('gibson.views',
	(r'^$', 'designer'),
	(r'^settings/$', 'construct_settings'),
	(r'^saveMeta/$', 'update_meta'),
	(r'^saveSettings/$', 'update_settings'),
	(r'^design/$', 'design_tab'),
)

urlpatterns = patterns('gibson.views',
	(r'^$', 'constructs'),
	(r'^add$', 'construct_add'),
	(r'^(?P<cid>\d+)/\w+\.gb', 'download'),
#	(r'^(?P<cid>\d+)/', include(constructpatterns)),
#	(r'^(?P<cid>\d+)/[\w\d ]+/', include(constructpatterns)),
	(r'^(?P<cid>\d+)/', include(designerpatterns)),
)

