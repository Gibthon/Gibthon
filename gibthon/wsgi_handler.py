import os, sys

sys.path.append('/home/bill/www')
sys.path.append('/home/bill/www/gibthon')

os.environ['DJANGO_SETTINGS_MODULE'] = 'gibthon.settings'

from django.core.handlers.wsgi import WSGIHandler
application = WSGIHandler()
