from django.conf.urls.defaults import patterns, include, url
from django.conf import settings
import views

urlpatterns = patterns('',
	(r'^$', views.apihelp),
	(r'^types/$', views.get_types),
	(r'^categories/', views.get_categories),
)
	