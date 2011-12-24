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
			'fragments[]': [f.id for f in c.fragments],
			'created': fragment.last_modified(),
		}
		return JsonResponse(ret)
	except ObjectDoesNotExist:
		return JsonResponse('Construct %s not found.' % s, ERROR)
	raise Http404

@login_required
def fragment_add(request, cid):
	con = get_construct(request.user, cid)
	if con:
		f = get_fragment(request.user, fid)
		if f:
			cf = con.add_fragment(f)
			if cf:
				if request.is_ajax():
					return JsonResponse({'fid': fid, 'cfid': cf.id,})
				return HttpResponseRedirect('/gibthon/%s/' % cid)
			else:
				if request.is_ajax():
					return JsonResponse('Could not add fragment %s to construct %s' % (fid, cid))
				raise Http404
		else:
			if request.is_ajax():
				return JsonResponse('Could not find fragment "%s"' % fid, ERROR)
			return HttpResponseNotFound()
	else:
		return HttpResponseNotFound()

@login_required
def fragment_delete(request, cid, cfid):
	con = get_construct(request.user, cid)
	if con:
		con.delete_cfragment(cfid)
		if request.is_ajax():
			return JsonResponse('OK');
		return HttpResponse("OK")
	else:
		if request.is_ajax():
			return JsonResponse('No such Construct ' + cid, ERROR)
		return HttpResponseNotFound()

@login_required
def save_order(request, cid):
	try:
		con = get_construct(request.user, cid)
		if request.method == 'POST':
			order = request.POST.getlist('order[]')
		if con and order:
			con.reorder_cfragments(order)
			return JsonResponse({'modified': con.last_modified(),});
		else:
			return HttpResponseNotFound()
	except Exception as e:
		print 'Error: %s' % (e.message)
	raise Http404
