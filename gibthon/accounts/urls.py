from django.conf.urls.defaults import patterns, include, url
from django.conf import settings

messagepatterns = patterns('accounts.views',
	(r'^detail$', 'message_detail'),
	(r'^delete$', 'message_delete'),
	(r'^add$', 'message_add'),
)

inboxpatterns = patterns('accounts.views',
	(r'^$', 'inbox'),
	(r'^all$', 'messages'),
	(r'^fetch$', 'fetch'),
	(r'^add$', 'message_add_all'),
	(r'^(?P<mid>\d+)/', include(messagepatterns)),
)

urlpatterns = patterns('accounts.views',
	(r'^$', 'redirect_home'),
	(r'^login/$', 'login'),
	(r'^logout$', 'logout'),
	(r'^profile/$', 'profile'),
#	(r'^reset/$', 'reset'),
	(r'^register/$', 'register'),
	(r'^register/(?P<email_hash>\w+)/$', 'create'),
	(r'^inbox/', include(inboxpatterns)),
)