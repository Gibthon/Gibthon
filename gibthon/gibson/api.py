import simplejson as json
from gibthon.jsonresponses import JsonResponse, ERROR

from fragment.models import *
from fragment.api import get_gene
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound, Http404
from gibson.views import get_construct
from fragment.views import get_fragment
from fragment.api import get_gene, read_meta

from django.core.exceptions import ObjectDoesNotExist
from django.contrib.auth.decorators import login_required

from forms import SettingsForm

def cf2dict(cf):
	"""
		Convert a ConstructFragment to a dictionary suitable for JSON encoding
	"""
	ret = {	
			'id': cf.id,
			'fid': cf.fragment.id,
			'order': cf.order,
			'direction': 1,
			's_feat': -1,
			's_offset': cf.start_offset,
			'e_feat': -1,
			'e_offset': cf.end_offset,
		}
	if cf.direction == 'r':
		ret['direction'] = -1

	if cf.start_feature:
		ret['s_feat'] = cf.start_feature.id
	if cf.end_feature:
		ret['e_feat'] = cf.end_feature.id
	
	return ret
			

@login_required
def save_meta(request, cid):
	"""	Save a construct's metadata
		Post Data, all optional
		name, desc
	"""
	if request.method == 'POST':
		con = get_construct(request.user, cid)
		if not con:
			return JsonResponse({'errors': 
				{'all': "Construct with id '%s' not found" % cid,}}, ERROR)
		name = request.POST.get('name', con.name)
		desc = request.POST.get('desc', con.description)
		print 'got name: "%s" desc: "%s"' %(name, desc)
		print 'post data: %s' % request.POST	
		if (name != con.name) or (desc != con.description):
			print 'saving name: "%s" desc: "%s"' %(name, desc)
			con.name = name
			con.description = desc
			try:
				con.save()
			except Exception as e:
				return JsonResponse('One or more fields are too long.', ERROR)
				
		return JsonResponse({'modified': con.last_modified(), 
			'fields': {'name': con.name, 'desc': con.description}});
		
	raise Http404
	
@login_required
def save_settings(request, cid):
	print 'update_settings request.method = %s' % request.method
	if request.method == 'POST':
		con = get_construct(request.user, cid)
		if not con:
			return JsonResponse({
				'errors': {
					'all':"Construct with id '%s' not found",
					},} % cid, ERROR)
		form = SettingsForm(request.POST, instance=con.settings)
		if form.is_valid():
			form.save()
			con.reset()
			data = {}
			for key,value in form.cleaned_data.items():
				data[key] = str(value);
			return JsonResponse({'modified': con.last_modified(), 'fields': data})
		return JsonResponse({'errors': form.errors,}, ERROR)
	raise Http404

@login_required
def get_info(request, cid):
	"""Get information about a construct"""
	try:
		c = get_construct(request.user, cid)
		cfs = []
		fs = []
		for cf in c.cf.all():
			cfs.append(cf2dict(cf));
			g = get_gene(request.user, cf.fragment.id)
			fs.append(read_meta(g))
			
		ret = {
            'id': cid,
			'name': c.name,
			'desc': c.description,
			'length': c.length(),
			'cfs': cfs,
			'fs': fs,
			'created': c.last_modified(),
		}
		return JsonResponse(ret)
	except ObjectDoesNotExist:
		return JsonResponse('Construct %s not found.' % s, ERROR)
	raise Http404

@login_required
def fragment_add(request, cid):
	con = get_construct(request.user, cid)
	if con:
		fid = request.POST.get('fid')
		
		if not fid:
			return JsonResponse("No fragment id provided", ERROR)
		f = get_fragment(request.user, fid)
		
		pos = request.POST.get('pos', 0);
		strand = request.POST.get('dir', 1);

		try:
			pos = int(pos)
		except ValueError as e:
			return JsonResponse('Position must be an integer, not "%s"' % pos, ERROR)
		try:
			strand = int(strand)
		except ValueError as e:
			return JsonResponse('Direction must be an integer, not "%s"' % strand,
					ERROR)
		if strand not in [-1, 1]:
			return JsonResponse('Strand must be 1 or -1, not "%s"' %strand, ERROR)
		
		direction = 'f'
		if strand == -1:
			direction = 'r'
		
		if f:
			cf = con.add_fragment(f, pos, direction)
			if cf:
				return JsonResponse(cf2dict(cf))
			else:
				return JsonResponse('Could not add fragment %s to construct %s' % (fid,
					cid), ERROR)
		else:
			if request.is_ajax():
				return JsonResponse('Could not find fragment "%s"' % fid, ERROR)
			return HttpResponseNotFound()
	else:
		return HttpResponseNotFound()

@login_required
def fragment_remove(request, cid):
	con = get_construct(request.user, cid)
	#post = json.loads(request.raw_post_data)
	cfid = request.POST.get('cfid')
	if not cfid:
		return JsonResponse('No ConstructFragment ID specified', ERROR)
	if con:
		try:
			con.delete_cfragment(cfid)
		except ObjectDoesNotExist as e:
			return JsonResponse('No such fragment "%s" associated with construct "%s".' % (cfid, cid))
		
		return JsonResponse('OK');
	else:
		return JsonResponse('No such Construct ' + cid, ERROR)

directions = {1:'f', -1:'r',}

@login_required
def save_order(request, cid):
	con = get_construct(request.user, cid)
	cfid = request.POST.getlist('cfid[]')
	direction = request.POST.getlist('direction[]')
	if not cfid:
		return JsonResponse('No cfids provided.', ERROR)
	if not direction:
		return JsonResponse('No directions provided', ERROR)
	if len(cfid) != len(con.cf.all()):
		return JsonResponse('Only %s ConstructFragments provided, %s required.' 
				% (len(cfid), len(con.cf.all())), ERROR)
	if len(direction) != len(con.cf.all()):
		return JsonResponse('Only %s directions provided, %s required.' 
				% (len(direction), len(con.cf.all())), ERROR)
	
	if con:
		dirs = []
		for d in direction:
			dirs.append(directions.get(d, ' '))
		con.reorder_cfragments(cfid, dirs)
		return JsonResponse({'modified': con.last_modified(),});
	else:
		return HttpResponseNotFound()
