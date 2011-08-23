from django.conf.urls.defaults import patterns, include, url
from django.conf import settings

messagepatterns = patterns('accounts.views',
	(r'^detail$', 'message_detail'),
)

urlpatterns = patterns('accounts.views',
	(r'^$', 'redirect_home'),
	(r'^login/$', 'login'),
	(r'^logout$', 'logout'),
	(r'^profile/$', 'profile'),
#	(r'^reset/$', 'reset'),
	(r'^register/$', 'register'),
	(r'^register/(?P<email_hash>\w+)/$', 'create'),
	(r'^inbox/$', 'inbox'),
	(r'^inbox/fetch$', 'fetch'),
	(r'^inbox/(?P<mid>\d+)/', include(messagepatterns)),
)