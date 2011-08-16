"""These views are simply a JSON API for the fragment view"""

from fragment.models import *
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound, Http404
from django.contrib.auth.decorators import login_required
from django.core.exceptions import ObjectDoesNotExist

import simplejson as json

OK = 0
ERROR = -1

class JsonResponse(HttpResponse):
	def __init__(self, data, state = OK):
		print "returning JSON: '%s'" % json.dumps([state, data]) 
		HttpResponse.__init__(self, json.dumps([state, data]), mimetype='application/json')

get_values = ['meta', 'seq', 'annotations', 'refs', 'feats']

@login_required
def get(request, _id):
	"""Handles a JS request for data"""
	try:
		the_id = int(_id)
	except ValueError:
		raise Http404
	print "api.get(request, %i) method: '%s', GET: '%s'" % (the_id, request.method, request.GET)
	#if request.is_ajax() and 'value' in request.GET:
	if 'value' in request.GET:
		value = request.GET['value'].lower()
		if not value in get_values:
			return JsonResponse("ERROR: Invalid value '%s'." % value, ERROR)
		try:
			g = Gene.objects.get(id = the_id, owner=request.user)
		except ObjectDoesNotExist:
			return JsonResponse("ERROR: Fragment with ID='%s' does not exist." % id, ERROR)
		
		if value == 'meta':
			data = {	'name': g.name,
						'desc': g.description,
						'origin': g.get_origin_display(),
					 }
		elif value == 'seq':
			#return a section of the sequence
			try:
				offset = int(request.GET.get('offset', 0))
			except ValueError:
				return JsonResponse("ERROR: Invalid offset '%s'." % request.GET.get('offset', 0), ERROR)
			try:
				length = int(request.GET.get('length', 1000))
			except ValueError:
				return JsonResponse("ERROR: Invalid length '%s'." % request.GET.get('length', 1000), ERROR)
			data = g.sequence[offset : offset+length]
		elif value == 'annotations':
			#return a dict of the annotations
			data = {}
			for a in g.annotations.all():
				if a.key not in data:
					data[a.key] = []
				data[a.key].append(a.value)
		elif value == 'refs':
			#get references
			data = []
			for r in g.references.all():
				data.append({	'title': r.title,
									'authors':r.authors,
									'journal':r.journal,
									'medline_id':r.medline_id,
									'pubmed_id':r.pubmed_id,
								})
		elif value == 'feats':
			#get features
			data = []
			for f in g.features.all():
				quals = []
				for q in f.qualifiers.all():
					quals.append({	'name': q.name,
										'data': q.data,
									 })
				s = None
				if f.direction == 'f':
					s = 1
				elif f.direction == 'r':
					s = -1
				data.append({	'start': f.start,
									'end': f.end,
									'strand': s,
									'type': f.type,
									'qualifiers': quals,
								})
		#all went well, return the JSON
		return JsonResponse(data)
	raise Http404