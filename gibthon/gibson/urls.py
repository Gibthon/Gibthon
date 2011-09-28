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

constructpatterns = patterns('gibson.designer',
	(r'^$', 'designer'),
	(r'^design/$', 'design_tab'),
	(r'^settings/$', 'construct_settings'),
	(r'^saveMeta/$', 'update_meta'),
	(r'^saveSettings/$', 'update_settings'),
)

urlpatterns = patterns('gibson.views',
	(r'^$', 'constructs'),
	(r'^add$', 'construct_add'),
	(r'^(?P<cid>\d+)/\w+\.gb', 'download'),
	(r'^(?P<cid>\d+)/delete/$', 'construct_delete'),
	(r'^(?P<cid>\d+)/fragments/$', 'construct_fragment'),
	(r'^(?P<cid>\d+)/addFragment/(?P<fid>\d+)/$', 'fragment_add'),
	(r'^(?P<cid>\d+)/', include(constructpatterns)),
)

