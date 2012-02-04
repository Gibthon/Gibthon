from django.conf.urls.defaults import patterns, include, url
from django.conf import settings

bufferpatterns = patterns( 'digest.views.buffers', 
	( r'^$', 'buffers' ),
	( r'group$', 'group' ),
	( r'manufacturer$', 'manufacturer' ),
)

enzymepatterns = patterns( 'digest.views.enzymes',
	( r'^$', 'enzymes' ),
	( r'manufacturers$', 'manufacturer' ),
)

urlpatterns = patterns('digest.views.base',
	( r'^$', 'base' ),
	( r'^buffers/', include( bufferpatterns ) ),
	( r'^enzymes/', include( enzymepatterns ) ),
)
