import os, sys


sys.path.append('/home/bill/www')
sys.path.append('/home/bill/www/gibthon')

os.environ['DJANGO_SETTINGS_MODULE'] = 'gibthon.settings'

import django.core.handlers.wsgi
application = django.core.handlers.wsgi.WSGIHandler()
