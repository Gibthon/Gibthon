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
)

apipatterns = patterns('gibson.api',
	(r'^saveSettings/$', 'update_settings'),
	(r'^saveMeta/$', 'save_meta'),
	(r'^getInfo/$', 'get_info'),
	(r'^addFragment/(?P<fid>\d+)/$', 'fragment_add'),
	(r'^removeFragment/(?P<cfid>\d+)/$', 'fragment_delete'),
	(r'^saveOrder/$', 'save_order'),
)

urlpatterns = patterns('gibson.views',
	(r'^$', 'constructs'),
	(r'^add$', 'construct_add'),
	(r'^(?P<cid>\d+)/\w+\.gb', 'download'),
	(r'^(?P<cid>\d+)/delete/$', 'construct_delete'),
	(r'^(?P<cid>\d+)/fragments/$', 'construct_fragment'),
	(r'^(?P<cid>\d+)/process/$', 'process'),
	(r'^(?P<cid>\d+)/clipping/(?P<cfid>\d+)/$', 'fragment_clipping'),
	(r'^(?P<cid>\d+)/clipping/(?P<cfid>\d+)/apply/$', 'apply_clipping'),
	(r'^(?P<cid>\d+)/summary/$', 'summary'),
	(r'^(?P<cid>\d+)/primers/', include(primerpatterns)),
	(r'^api/(?P<cid>\d+)/', include(apipatterns)),
	(r'^(?P<cid>\d+)/', include(constructpatterns)),
)

