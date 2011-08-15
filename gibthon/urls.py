from django.conf.urls.defaults import patterns, include, url
from django.conf import settings
from django.contrib import admin
admin.autodiscover()


userpatterns = patterns('gibthon.views',
	(r'^$', 'redirect_home'),
	(r'^login/$', 'login'),
	(r'^logout$', 'logout'),
	(r'^profile/$', 'profile'),
#	(r'^reset/$', 'reset'),
	(r'^register/$', 'register'),
	(r'^register/(?P<email_hash>\w+)/$', 'create'),
)

toolpatterns = patterns('gibthon.views',
	(r'^$', 'redirect_home'),
	(r'^ligate/$', 'ligcal'),
	(r'^ligcalc/$', 'ligate'),
	(r'^molcalc/$', include('gibthon.molcal.urls')),
	(r'^molcal/$', 'molcal'),
)

urlpatterns = patterns('',
	(r'^$', 'gibthon.views.home'),
	(r'^gibthon/', include('gibthon.gibson.urls')),
	(r'^fragment/', include('gibthon.fragment.urls')),
	(r'^help/', include('gibthon.help.urls')),
    (r'^user/', include(userpatterns)),
    (r'^tools/', include(toolpatterns)),
    (r'^admin/', include(admin.site.urls)),
    (r'^admin/doc/', include('django.contrib.admindocs.urls')),
    (r'^captcha/', include('captcha.urls')),
    (r'^partsearch/', include('partsearch.urls')),
)

