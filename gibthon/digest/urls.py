from django.conf.urls.defaults import patterns, include, url
from django.conf import settings

bufferpatterns = patterns( 'digest.views', 
	( r'^$', 'buffers' ),
	( r'group$', 'group' ),
	( r'manufacturer$', 'manufacturer' ),
)

enzymepatterns = patterns( 'digest.views',
	( r'^$', 'enzymes' ),
)

urlpatterns = patterns('digest.views',
	( r'^$', 'base' ),
	( r'^buffers/', include( bufferpatterns ) ),
	( r'^enzymes/', include( enzymepatterns ) ),
)
