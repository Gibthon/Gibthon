
############### JSON API for designer

from gibthon.jsonresponses import JsonResponse

@login_required
def update_meta(request, cid): # 'saveMeta/'
	if request.method == 'POST': #and request.is_ajax():
		con = get_construct(request.user, cid)
		if not con:
			return JsonResponse({'errors': {'all': "Construct with id '%s' not found" % cid,}}, ERROR)
		name = request.POST.get('name', con.name)
		desc = request.POST.get('desc', con.description)
		if (name != con.name) or (desc != con.description):
			con.name = name
			con.description = desc
			con.save()
		
		return JsonResponse({'modified': con.last_modified(), 'fields': {'name': name, 'desc': desc}});
		
	raise Http404
	
@login_required
def update_settings(request, cid):
	if request.method == 'POST': #and request.is_ajax():
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
	

######################################### end JSON
