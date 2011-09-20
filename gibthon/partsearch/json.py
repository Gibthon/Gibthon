from django.http import HttpResponse

import simplejson as json

class JsonpResponse(HttpResponse):
	def __init__(self, data, request):
		callback = request.GET.get('callback', 'callback')
		json_data = json.dumps(data)
		jsonp = "%s(%s)" % (callback, json_data)
		HttpResponse.__init__(self, jsonp, mimetype='application/json')
