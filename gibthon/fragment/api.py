"""These views are simply a JSON API for the fragment view"""

from fragment.models import *
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound, Http404
from django.shortcuts import render_to_response
from django.contrib.auth.decorators import login_required
from django.core.exceptions import ObjectDoesNotExist

import simplejson as json

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

OK = 0
ERROR = -1

class JsonResponse(HttpResponse):
	def __init__(self, data, state = OK):
		HttpResponse.__init__(self, json.dumps([state, data]), mimetype='application/json')

class RawJsonResponse(HttpResponse):
	def __init__(self, data):
		print "Return JSON: '%s'" % json.dumps(data)
		HttpResponse.__init__(self, json.dumps(data), mimetype='application/json')

# functions which get the appropriate data
def get_meta(g, request):
	return JsonResponse({	'name': g.name,
									'desc': g.description,
									'origin': g.get_origin_display(),
									'length': len(g.sequence)
								})

def get_seq_meta(g, request):
	"""get all the sequence metadata"""		
	#get features
	feats = []
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
		feats.append({	'start': f.start,
							'end': f.end,
							'strand': s,
							'type': f.type,
							'qualifiers': quals,
						})
		
	#assume Ambiguous DNA
	let = Seq(IUPAC.IUPACAmbiguousDNA.letters, IUPAC.IUPACAmbiguousDNA())
	rlet = let.complement()
	alpha = {}
	for i in range(len(let)):
		alpha[let[i].lower()] = rlet[i].lower()
		alpha[let[i].upper()] = rlet[i].upper()
	
	return JsonResponse({	'len': len(g.sequence),
									'feats': feats,
									'alpha': alpha,
								})
	
def get_seq(g, request):
	#return a section of the sequence
	try:
		offset = int(request.GET.get('offset', 0))
	except ValueError:
		return JsonResponse("ERROR: Invalid offset '%s'." % request.GET.get('offset', 0), ERROR)
	try:
		length = int(request.GET.get('length', 1000))
	except ValueError:
		return JsonResponse("ERROR: Invalid length '%s'." % request.GET.get('length', 1000), ERROR)
	return JsonResponse(g.sequence[offset : offset+length])
	
def get_annotations(g, request):
	#return a dict of the annotations
	data = {}
	for a in g.annotations.all():
		if a.key not in data:
			data[a.key] = []
		data[a.key].append(a.value)
	return JsonResponse(data)
	
def get_refs(g, request):
	#get references
	data = []
	for r in g.references.all():
		data.append({	'title': r.title,
							'authors':r.authors,
							'journal':r.journal,
							'medline_id':r.medline_id,
							'pubmed_id':r.pubmed_id,
						})
	return JsonResponse(data)
	
def get_feats(g, request):
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
	return JsonResponse(data)
	
def get_len(g,request):
	return JsonResponse(len(g.sequence))

def get_alpha(g, request):
	#assume Ambiguous DNA
	let = Seq(IUPAC.IUPACAmbiguousDNA.letters, IUPAC.IUPACAmbiguousDNA())
	rlet = let.complement()
	data = {}
	for i in range(len(let)):
		data[let[i].lower()] = rlet[i].lower()
		data[let[i].upper()] = rlet[i].upper()
	return JsonResponse(data)

get_map = 	{	'meta': get_meta,
				'seq': get_seq,
				'annotations': get_annotations,
				'refs': get_refs,
				'feats': get_feats,
				'len': get_len,
				'alpha': get_alpha,
				'seq_meta': get_seq_meta,
			}

@login_required
def get(request, _id):
	"""Handles a JS request for data"""
	try:
		the_id = int(_id)
	except ValueError:
		raise Http404

	#if request.is_ajax() and 'value' in request.GET:
	if 'value' in request.GET:
		value = request.GET['value'].lower()
		if not value in get_map:
			return JsonResponse("ERROR: Invalid value '%s'." % value, ERROR)
		try:
			g = Gene.objects.get(id = the_id, owner=request.user)
		except ObjectDoesNotExist:
			return JsonResponse("ERROR: Fragment with ID='%s' does not exist." % id, ERROR)
		return get_map[value](g, request)
	raise Http404


