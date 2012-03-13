"""Server-side Fragment Library API"""
#get/set = view
#read/write = database operation

from fragment.models import *
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound, Http404
from django.shortcuts import render_to_response
from django.contrib.auth.decorators import login_required
from django.core.exceptions import ObjectDoesNotExist

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

import simplejson as json
from gibthon.jsonresponses import JsonResponse, RawJsonResponse, ERROR

## Helpful functions

the_alphabet = {
	'A': 'T',
	'B': 'V',
	'C': 'G',
	'D': 'H',
	'G': 'C',
	'H': 'D',
	'K': 'M',
	'M': 'K',
	'N': 'N',
	'R': 'Y',
	'S': 'S',
	'T': 'A',
	'V': 'B',
	'W': 'W',
	'Y': 'R',
	'a': 't',
	'b': 'v',
	'c': 'g',
	'd': 'h',
	'g': 'c',
	'h': 'd',
	'k': 'm',
	'm': 'k',
	'n': 'n',
	'r': 'y',
	's': 's',
	't': 'a',
	'v': 'b',
	'w': 'w',
	'y': 'r',
}

def get_alphabet():
	return the_alphabet

def get_gene(usr, fid):
	try:
		fid = int(fid)
	except ValueError:
		raise Http404
	return Gene.objects.get(id = fid, owner=usr)
	
def read_min_meta(g):
	"""Return JSON-ready minimal metadata about a fragment"""
	f = {
			'id': g.id,
			'name': g.name,
			'desc': g.description,
			'origin': g.origin,
			'length': len(g.sequence),
		}

def read_meta(g):
	"""Return JSON-ready metadata about a fragment"""
	refs = []
	for r in g.references.all():
		refs.append( {
			'title': r.title,
			'authors': r.authors,
			'journal': r.journal,
			'medline_id': r.medline_id,
			'pubmed_id': r.pubmed_id,
		})
	annots = {}
	for a in g.annotations.all():
		annots[a.key] = a.value
	feats = []	
	for f in g.features.all():
		d = 1
		if f.direction == 'r':
			d = -1
		feat =  {
			'id': f.id,
			'type': f.type,
			'start': f.start,
			'end': f.end,
			'strand': d,
			'qualifiers': [],			
		}
		for q in f.qualifiers.all():
			feat['qualifiers'].append({'name': q.name, 'value':q.data,})
		feats.append(feat)

	return {	'name': g.name,
				'fid': g.id,
				'desc': g.description,
				'refs': refs,
				'annots': annots,
				'origin': g.get_origin_display(),
				'length': len(g.sequence),
				'feats': feats,
	}
	
def write_meta(g, meta):
	"""save meta to g"""
	g.name = meta.get('name', g.name)
	g.description = meta.get('desc', g.description)
	#only change references if they've been given
	if meta.has_key('refs'):
		#clear out previous references
		Reference.remove(g)
		Reference.add(g, meta['refs'])
	if meta.has_key('annots'):
		#clear out previous annotations
		Annotation.remove(g)
		for (key, value) in meta['annots']:
			Annotation.add(g, key, value)
	try:
		g.save()
	except Exception as e:
		print "raised exception of type %s: %s" % (type(e), e)
	
def read_features(g):
	"""Read all the features of a fragment"""
	data = []
	for f in g.features.all():
		quals = []
		for q in f.qualifiers.all():
			quals.append({	
				'name': q.name,
				'data': q.data,
			})
		s = None
		if f.direction == 'f':
			s = 1
		elif f.direction == 'r':
			s = -1
		data.append({
			'start': f.start,
			'end': f.end,
			'strand': s,
			'type': f.type,
			'id': f.id,
			'qualifiers': quals,
		})
	return data

def write_features(g, feats):
	"""Update the features in the database"""
	#remove old features
	Feature.remove(g)
	for f in feats:
		#make sure data makes sense
		start = int(f.get('start',0))
		end = int(f.get('end',0))
		strand = int(f.get('strand', 1))
		if strand not in [-1, 1]:
			strand = 1
		if strand == -1 and start < end:
			t = end
			end = start
			start = t
		elif strand == 1 and start > end:
			t = end
			end = start
			start = t
		
		quals = {}
		for q in f.get('qualifiers'):
			quals[q.get('name')] = q.get('value')
		
		ft = SeqFeature(
			location=FeatureLoaction(start, end), 
			type=f.get('type', ''), 
			strand=strand,
			id=f.get('id'), 
			qualifiers=quals,
		)
			
		Feature.add(ft, g)
	
# Actual API stuff

@login_required
def list_fragments(request):
	"""Return metadata of all fragments owned by a user"""
	try:
		frags = Gene.objects.filter(owner=request.user)
	except ObjectDoesNotExist:
		frags = []
	ret = []
	for f in frags:
		ret.append(read_meta(f))
	return JsonResponse(ret)

@login_required
def get_meta(request, fid):
	"""get a particular fragment's metadata"""
	g = get_gene(request.user, fid)
	return JsonResponse( read_meta(g))

@login_required
def get_multi_meta(request):
	"""return metadata for several fragments"""
	fids = request.POST.getlist('fids[]')
	if fids:
		r = {}
		for fid in fids:
			g = get_gene(request.user, fid)
			if g:
				r[fid] = read_meta(g)
		return JsonResponse(r)
	raise Http404

@login_required
def set_meta(request, fid):
	"""Update a fragment's metadata"""
  
	if request.method == 'POST':
		try:
			g = get_gene(request.user, fid)
			meta = json.loads(request.raw_post_data)
			if meta:
				write_meta(g, meta)
				return get_meta(request, fid)
			return JsonResponse("No metadata supplied.", ERROR)			
		except ObjectDoesNotExist:
			return JsonResponse("Fragment with ID='%s' does not exist." % fid, ERROR)
	raise Http404

@login_required
def get_features(request, fid):
	"""Get a fragment's features"""
	g = get_gene(request.user, fid)
	return JsonResponse({
		'feats': read_features(g),
		'alpha': get_alphabet(),
		'length': len(g.sequence),
	})

@login_required
def set_features(request, fid):
	"""Save a fragment's features"""
	if request.method == 'POST':
		try:
			g = get_gene(request.user, fid)
			feats = request.POST.get('features')
			if feats:
				write_features(g, feats)
				return JsonResponse('Done')
			return JsonResponse('No features provided', ERROR)
		except ObjectDoesNotExist:
			JsonResponse("Fragment with ID='%s' does not exist" % fid, ERROR)
	raise Http404

@login_required
def get_length(request, fid):
	g = get_gene(request.user, fid)
	return JsonResponse(len(g.seq))

@login_required
def get_sequence(request, fid):
	"""return a section of the sequence"""
	try:
		offset = int(request.GET.get('offset', 0))
	except ValueError:
		return JsonResponse("ERROR: Invalid offset '%s'." % request.GET.get('offset', 0), ERROR)
	try:
		length = int(request.GET.get('length', 1000))
	except ValueError:
		return JsonResponse("ERROR: Invalid length '%s'." % request.GET.get('length', 1000), ERROR)
	g = get_gene(request.user, fid)
	if length <= 0:
		return JsonResponse({
			'sequence':g.sequence[offset:],
			'alpha': get_alphabet(),
		})
	return JsonResponse(g.sequence[offset : offset+length])

# Editing of sequence data not yet supported
#def set_sequence(request, fid):
#	"""set the sequence"""
#	try:
#		offset = int(request.GET.get('offset', 0))
#	except ValueError:
#		return JsonResponse("ERROR: Invalid offset '%s'." % request.GET.get('offset', 0), ERROR)
	
