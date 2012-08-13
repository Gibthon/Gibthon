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
			'title':'Construct Designer',
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
			'cid': cid,
			'settings':s,
		})
		return HttpResponse(t.render(c))
	else:
		return HttpResponseNotFound()
