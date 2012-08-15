from django.conf.urls.defaults import patterns, include, url
from django.conf import settings
from django.contrib.auth.views import (password_reset, password_reset_done,
	password_reset_confirm, password_reset_complete)

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

passwordpatterns = patterns('django.contrib.auth.views',
	(r'^reset/$', 'password_reset', {
		'template_name': 'user/password_reset_form.html', 
		'email_template_name':'user/password_reset_email.html',
		'from_email': 'admin@gibthon.org',
		'post_reset_redirect' : '/user/password/reset/done/'
	}),
	(r'^reset/done/$', 'password_reset_done', {
		'template_name': 'user/password_reset_done.html'
	}),
	(r'^reset/(?P<uidb36>[0-9A-Za-z]+)-(?P<token>.+)/$',
		'password_reset_confirm', {
			'template_name': 'user/password_reset_confirm.html',
			'post_reset_redirect' : '/user/password/done/'
	}),
	(r'^done/$',	'password_reset_complete', {
		'template_name': 'user/password_reset_complete.html',
	}),
)

urlpatterns = patterns('accounts.views',
	(r'^$', 'redirect_home'),
	(r'^login/$', 'login'),
	(r'^logout$', 'logout'),
	(r'^profile/$', 'profile'),
#	(r'^reset/$', 'reset'),
	(r'^register/$', 'register'),
	(r'^register/(?P<email_hash>\w+)/$', 'create'),
	(r'^delete/$', 'delete'),
	(r'^inbox/', include(inboxpatterns)),
	(r'^password/', include(passwordpatterns)),
)
