import simplejson as json
from gibthon.jsonresponses import JsonResponse, ERROR

from gibson.views import get_construct

from django.core.exceptions import ObjectDoesNotExist

@login_required
def save_meta(request, cid):
	"""	Save a construct's metadata
		Post Data, all optional
		name, desc
	"""
	if request.method == 'POST':
		con = get_construct(request.user, cid)
		if not con:
			return JsonResponse({'errors': {'all': "Construct with id '%s' not found" % cid,}}, ERROR)
		name = request.POST.get('name', con.name)
		desc = request.POST.get('desc', con.description)
		if (name != con.name) or (desc != con.description):
			con.name = name
			con.description = desc
			try:
				con.save()
			except Exception as e:
				return JsonResponse('One or more fields are too long.', ERROR)
				
		return JsonResponse({'modified': con.last_modified(), 'fields': {'name': con.name, 'desc': con.desc}});
		
	raise Http404
	
@login_required
def update_settings(request, cid):
	if request.method == 'POST':
		con = get_construct(request.user, cid)
		if not con:
			return JsonResponse({'errors': {'all':"Construct with id '%s' not found",},} % cid, ERROR)
		form = SettingsForm(request.POST, instance=con.settings)
		if form.is_valid():
			form.save()
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
		c = get_construct(cid)
		ret = {
			'name': c.name,
			'desc': c.description,
			'length': c.length(),
			'fragments[]': [cf.fragment for cf in c.cf.all()],
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
		if f:
			cf = con.add_fragment(f)
			if cf:
				return JsonResponse({'fid': fid, 'cfid': cf.id,})
			else:
				return JsonResponse('Could not add fragment %s to construct %s' % (fid, cid))
		else:
			if request.is_ajax():
				return JsonResponse('Could not find fragment "%s"' % fid, ERROR)
			return HttpResponseNotFound()
	else:
		return HttpResponseNotFound()

@login_required
def fragment_remove(request, cid):
	con = get_construct(request.user, cid)
	fid = request.POST.get('fid')
	if not fid:
		return JsonResponse('No fragment ID specified', ERROR)
	if con:
		try:
			f = con.cf.get(fragment=fid)
			con.delete_cfragment(f.cf.id)
		except ObjectDoesNotExist as e:
			return JsonResponse('No such fragment "%s" associated with construct "%s".' % (fid, cid))
		
		return JsonResponse('OK');
	else:
		return JsonResponse('No such Construct ' + cid, ERROR)

directions = {1:'f', -1:'r',}

@login_required
def save_order(request, cid):
	try:
		con = get_construct(request.user, cid)
		order = request.POST.getlist('order[]')
		direction = request.POST.getlist('direction[]')
		if not order:
			return JsonResponse('No order provided.', ERROR)
		if len(order) != len(con.cf.all()):
			return JsonResponse('Only %s fragment ids provided, %s required.' % (len(order), len(con.cf.all())), ERROR)
		if not direction:
			direction = [0]*len(order)
		if len(direction) != len(order):
			return JsonResponse('len(direction) = %s != len(order) = %s' % (len(direction), len(order)))
		
		if con:
			t_dir = []
			for d in direction:
				if d in directions.keys():
					t_dir.append(directions[d])
				else:
					t_dir.append(' ')
			con.reorder_cfragments(order, t_dir)
			return JsonResponse({'modified': con.last_modified(),});
		else:
			return HttpResponseNotFound()
	except Exception as e:
		print 'Error: %s' % (e.message)
	raise Http404
