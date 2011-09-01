import os, sys

from localsettings import localsettings

sys.path.append(localsettings.gibthon_path)
sys.path.append(localsettings.gibthon_path + '/gibthon')

os.environ['DJANGO_SETTINGS_MODULE'] = 'gibthon.settings'

from django.core.handlers.wsgi import WSGIHandler
application = WSGIHandler()
