"""
Classes for JSON responses
"""

from django.http import HttpResponse
import simplejson as json

OK = 0
ERROR = -1

class JsonResponse(HttpResponse):
	def __init__(self, data, state = OK):
		HttpResponse.__init__(self, json.dumps([state, data]), mimetype='application/json')

class RawJsonResponse(HttpResponse):
	def __init__(self, data):
		HttpResponse.__init__(self, json.dumps(data), mimetype='application/json')

