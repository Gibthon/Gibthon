from partsearch.models import *
from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound, Http404
from django.core.exceptions import *
from django.shortcuts import render_to_response

from json import JsonpResponse

def apihelp(request):
	t = loader.get_template('partsearch/apihelp.html')
	c = RequestContext(request, {'title': 'PartSearch API',})
	return HttpResponse(t.render(c))

def get_types(request):
	types = PartType.objects.all();
	ret = []
	for t in types:
		ret.append(t.name)
	return JsonpResponse(ret, request)

def get_categories(request):
	cats = Category.objects.all();
	type = request.GET.get('type', 'array').lower()
	
	if type == 'array':
		ret = []
		for c in cats:
			ret.append(c.name)
	elif type == 'map':
		def add(d, dirs):	
			if dirs[0] not in d:
				d[dirs[0]] = {}
			if len(dirs) != 1:
				add(d[dirs[0]], dirs[1:])
		ret = {}
		for c in cats:
			add(ret, c.as_list())
				
	return JsonpResponse(ret, request)
