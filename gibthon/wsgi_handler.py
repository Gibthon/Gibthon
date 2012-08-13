import os, sys

import site
site.addsitedir('/usr/local/pythonenv/gibthon-env/lib/python2.7/site-packages')
site.addsitedir('/usr/local/pythonenv/gibthon-env/lib64/python2.7/site-packages')

gibthon_path = '/home/gibthon/site'

sys.path.append(gibthon_path)
sys.path.append(gibthon_path + '/gibthon')

os.environ['DJANGO_SETTINGS_MODULE'] = 'gibthon.settings'

from django.core.handlers.wsgi import WSGIHandler
application = WSGIHandler()
