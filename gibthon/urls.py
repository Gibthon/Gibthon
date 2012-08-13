from django.conf.urls.defaults import patterns, include, url
from django.conf import settings
from django.contrib import admin
from django.views.generic.simple import direct_to_template, redirect_to
from django.contrib.staticfiles.urls import staticfiles_urlpatterns


admin.autodiscover()


toolpatterns = patterns('gibthon.views',
	(r'^$', 'redirect_home'),
	(r'^ligate/$', 'ligcal'),
	(r'^ligcalc/', include('gibthon.ligcalc.urls')),
	(r'^molcalc/', include('gibthon.molcal.urls')),
	(r'^molcal/$', 'molcal'),
	(r'^digest/', include( 'gibthon.digest.urls' ) ),
)

urlpatterns = patterns('',
	(r'^$', direct_to_template, {'template':'home.html'}),
	(r'^changes/$', 'gibthon.views.changes'),
	(r'^about/$', 'gibthon.views.about'),
	(r'^gibthon/', include('gibthon.gibson.urls')),
	(r'^fragment/', include('gibthon.fragment.urls')),
	(r'^help/', include('gibthon.help.urls')),
  (r'^user/', include('gibthon.accounts.urls')),
  (r'^tools/', include(toolpatterns)),
  (r'^admin/', include(admin.site.urls)),
  (r'^admin/doc/', include('django.contrib.admindocs.urls')),
  (r'^captcha/', include('captcha.urls')),
  (r'^partsearch/', include('partsearch.urls')),
	(r'pricalc/', include('gibthon.pricalc.urls')),
  (r'^robots.txt$', direct_to_template, {
		'template': 'robots.txt',
		'mimetype':'text/plain'
	}),
  (r'^sitemap.xml$', direct_to_template, {
		'template': 'sitemap.xml',
		'mimetype':'text/xml'
	}),
)


if settings.DEBUG:
	urlpatterns += patterns('',
		url(r'^media/(?P<path>.*)$', 'django.views.static.serve', {
		'document_root': settings.MEDIA_ROOT,
		}),
	)
	urlpatterns += staticfiles_urlpatterns()


