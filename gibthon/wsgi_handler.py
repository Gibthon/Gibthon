import os, sys

gibthon_path = '/home/bill/www'

sys.path.append(gibthon_path)
sys.path.append(gibthon_path + '/gibthon')

os.environ['DJANGO_SETTINGS_MODULE'] = 'gibthon.settings'

from django.core.handlers.wsgi import WSGIHandler
application = WSGIHandler()
