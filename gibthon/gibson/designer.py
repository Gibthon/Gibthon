from django.contrib.auth.decorators import login_required
from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound, Http404
from views import get_construct

from models import *
from forms import *

@login_required
def designer(request,cid):
	con = get_construct(request.user, cid)
	if con:
		t = loader.get_template('gibson/designer.html')
		c = RequestContext(request,{
			'title':'Construct Designer: '+con.name+'',
			'id': cid,
			'construct':con,
		})
		return HttpResponse(t.render(c))
	else:
		raise Http404()

@login_required
def design_tab(request, cid):
	con = get_construct(request.user, cid)
	if con:
		t = loader.get_template('gibson/designtab.html')
		c = RequestContext(request,{
			'id': cid,
			'construct':con,
		})
		return HttpResponse(t.render(c))
	else:
		return HttpResponseNotFound()

@login_required
def construct_settings(request, cid):
	con = get_construct(request.user, cid)
	if con:
		if request.method == 'POST':
			form = SettingsForm(request.POST, instance=con.settings)
			if form.is_valid():
				form.save()
				return HttpResponse()
		t = loader.get_template('gibson/settings.html')
		s = con.settings
		c = RequestContext(request, {
			'settings':s,
		})
		return HttpResponse(t.render(c))
	else:
		return HttpResponseNotFound()


############### JSON API for designer

from gibthon.jsonresponses import JsonResponse, ERROR

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
