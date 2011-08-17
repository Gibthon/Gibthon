from models import *
from forms import *
from fragment.models import Gene, Feature
from fragment.views import get_fragment

from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound
from django.core.context_processors import csrf
from django.contrib.auth.decorators import login_required
from django.views.decorators.http import condition
from django.core.exceptions import *

import csv
import time
import json
from copy import copy

def fix_request(reqp):
	rp = copy(reqp)
	for k,v in rp.iteritems():
		if k.startswith('shape'):
			rp['shape'] = v
			del rp[k]
			return rp
		else:
			continue

def get_construct(user, cid):
	try:
		c = Construct.objects.get(pk=cid, owner=user)
	except ObjectDoesNotExist:
		return False
	else:
		return c

def get_primer(user, pid):
	try:
		p = Primer.objects.get(pk=pid)
	except ObjectDoesNotExist:
		return False
	else:
		return p

@login_required
def download(request, cid):
	con = get_construct(request.user, cid)
	if con:
		#response = HttpResponse(mimetype='chemical/seq-na-genbank')
		response = HttpResponse(mimetype='text/plain')
		response.write(con.gb())
		return response
	else:
		return HttpResponseNotFound()

@login_required
def settings(request, cid):
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

@login_required
def settings_edit(request, cid):
	con = get_construct(request.user, cid)
	if con:
		t = loader.get_template('gibson/settings_edit.html')
		f = SettingsForm(instance=con.settings)
		f.is_valid()
		c = RequestContext(request, {
			'construct':con,
			'form':f,
		})
		return HttpResponse(t.render(c))
	else:
		return HttpResponseNotFound()
		
@login_required
def primer(request, cid, pid):
	con = get_construct(request.user, cid)
	if con:
		t = loader.get_template('gibson/primers.html')
		p = get_primer(request.user, pid)
		if p:
			c = RequestContext(request, {
				'construct': con,
				'primer_list': con.primer.all(),
				'title': 'Primers('+con.name+')',
				'primer': p.id,
			})
			return HttpResponse(t.render(c))
		else:
			return HttpResponseNotFound()
	else:
		return HttpResponseNotFound()

@login_required
def csv_primers(request, cid):
	con = get_construct(request.user, cid)
	if con:
		resp = HttpResponse(mimetype="text/csv")
		resp['Content-Disposition'] = 'attachment; filename=' + con.name + '-primers.csv'
		writer = csv.writer(resp)
		writer.writerow(['Name', 'Length', 'Melting Temperature', 'Sequence'])
		for p in con.primer.all():
			writer.writerow(p.csv())
		return resp
	else:
		return HttpResponseNotFound()

@login_required
def load_primer(request, cid, pid):
	con = get_construct(request.user, cid)
	if con:
		t = loader.get_template('gibson/primer.html')
		p = get_primer(request.user, pid)
		if p:
			c = RequestContext(request, {
				'primer': p,
			})
			return HttpResponse(t.render(c))
		else:
			return HttpResponseNotFound()
	else:
		return HttpResponseNotFound()
		
@login_required
def primers(request, cid):
	con = get_construct(request.user, cid)
	if con and len(con.primer.all()) == 2*len(con.cf.all()):
		t = loader.get_template('gibson/primers.html')
		c = RequestContext(request, {
			'construct': con,
			'primer_list': con.primer.all(),
			'title':'Primers('+con.name+')',
		})
		return HttpResponse(t.render(c))
	else:
		return HttpResponseNotFound()

@login_required
@condition(etag_func=None)
def process(request, cid, reset=True, new=True):
	con = get_construct(request.user, cid)
	if con:
		resp = HttpResponse(con.process(reset, new), mimetype="text/plain")
		return resp
	else:
		return HttpResponseNotFound()

@login_required
def constructs(request):
	t = loader.get_template('gibson/constructs.html')
	con = Construct.objects.all().filter(owner=request.user)
	form = ConstructForm()
	c = RequestContext(request, {
		'construct_list':con,
		'title':'Construct Designer',
		'construct_form':form,
	})
	c.update(csrf(request))
	return HttpResponse(t.render(c))

@login_required	
def construct_add(request):
	if request.method == 'POST':
		rp = fix_request(request.POST)
		form = ConstructForm(rp)
		if form.is_valid():
			c = form.save()
			c.owner = request.user
			c.save()
			return HttpResponse('/gibthon/'+str(c.id)+'/'+c.name+'/')
		t = loader.get_template('gibson/gibsonnew.html')
		con = Construct.objects.all().filter(owner=request.user)
		c = RequestContext(request, {
			'construct_list':con,
			'title':'Construct Designer',
			'construct_form':form,
		})
		c.update(csrf(request))
		return HttpResponse(t.render(c))
	else:
		return HttpResponseNotFound()

@login_required
def construct(request,cid):
	con = get_construct(request.user, cid)
	if con:
		t = loader.get_template('gibson/gibson.html')
		if request.method == 'POST':
			rp = fix_request(request.POST)
			f = ConstructForm(rp, instance=con)
			f.save()
			return HttpResponse('/gibthon/'+str(con.id)+'/'+con.name+'/')
		construct_form = ConstructForm(instance=con)
		fragment_list = con.fragments.all().order_by('cf')
		feature_list = [FeatureListForm(f, con) for f in fragment_list]
		cf_list = con.cf.all()
		list = zip(fragment_list, feature_list, cf_list)
		c = RequestContext(request,{
			'title':'Construct Designer('+con.name+')',
			'list':list,
			'construct':con,
			'construct_form':construct_form,
		})
		return HttpResponse(t.render(c))
	else:
		return HttpResponseNotFound()

@login_required
def construct_delete(request, cid):
	con = get_construct(request.user, cid)
	if con:
		con.delete()
		return HttpResponseRedirect('/gibthon')
	else:
		return HttpResponseNotFound()

@login_required	
def fragment_viewer(request, cid, fid):
	if get_construct(request.user, cid):
		f = get_fragment(request.user, fid)
		if f:
			t = loader.get_template('gibson/fragment_viewer.html')
			c = RequestContext(request,{
				'fragment':f,
			})
			return HttpResponse(t.render(c))
		else:
			return HttpResponseNotFound()
	else:
		return HttpResponseNotFound()
	
@login_required
def fragment_browse(request, cid):
	con = get_construct(request.user, cid)
	if con:
		ex = [e.pk for e in con.fragments.all()]
		f = Gene.objects.all().exclude(pk__in=ex).filter(owner=request.user)
		t = loader.get_template('gibson/fragment_browser.html')
		c = RequestContext(request,{
			'fragment_list':f,
			'construct_id':cid,
		})
		return HttpResponse(t.render(c))
	else:
		return HttpResponseNotFound()
	
@login_required
def fragment_add(request, cid, fid):
	con = get_construct(request.user, cid)
	if con:
		f = get_fragment(request.user, fid)
		if f:
			add_fragment(con,f)
			return HttpResponseRedirect('/gibthon/'+cid+'/'+con.name+'/')
		else:
			return HttpResponseNotFound()
	else:
		return HttpResponseNotFound()

@login_required
def fragment_delete(request, cid, fid):
	con = get_construct(request.user, cid)
	if con:
		try: cf = ConstructFragment.objects.get(fragment=fid, construct=cid)
		except ObjectDoesNotExist: return HttpResponseNotFound()
		else:
			cf.delete()
			return HttpResponseRedirect('/gibthon/'+cid+'/'+cf.construct.name+'/')
	else:
		return HttpResponseNotFound()

@login_required
def save(request, cid):
	con = get_construct(request.user, cid)
	if request.method == 'POST' and con:
		t = loader.get_template('gibson/date.html')
		order = request.POST.getlist('order[]')
		feature_select = request.POST.getlist('feature_select[]')
		direction = request.POST.getlist('direction[]')
		for i,cff in enumerate(zip(order,feature_select,direction)):
			[cfid,fids, d] = cff
			cf = ConstructFragment.objects.get(pk=cfid)
			cf.order = i
			[fsid, feid] = fids.split(',')
			cf.start_feature = Feature.objects.get(pk=fsid)
			cf.end_feature = Feature.objects.get(pk=feid)
			cf.direction = d
			cf.save()
		con.save()
		c = RequestContext(request,{'date':con.modified})
		return HttpResponse(t.render(c))
	else:
		return HttpResponseNotFound()
		
@login_required
def summary(request, cid):
	con = get_construct(request.user, cid)
	t = loader.get_template('gibson/summary.html')
	c = RequestContext(request, {
		'features':con.features_pretty(),
	})
	return HttpResponse(t.render(c))

@login_required
@condition(etag_func=None)
def primer_offset(request, cid, pid):
	con = get_construct(request.user, cid)
	p = get_primer(request.user, pid)
	if request.method == 'POST' and con and p:
		if 'flap' in request.POST:
			offset = request.POST['flap']
			if p.flap.top:
				p.flap.cfragment.start_offset = offset
			else:
				p.flap.cfragment.end_offset = offset
			p.flap.cfragment.save()
		else:
			offset = request.POST['stick']
			if p.flap.top:
				p.stick.cfragment.end_offset = offset
			else:
				p.stick.cfragment.start_offset = offset
			p.stick.cfragment.save()
		return process(request, cid, reset=False, new=False)		
	else:
		return HttpResponseNotFound()

@login_required
@condition(etag_func=None)
def primer_reset(request, cid):
	con = get_construct(request.user, cid)
	if request.method == 'GET' and con:
		return process(request, cid, reset=True, new=False)
	else:
		return HttpResponseNotFound()

@login_required
def primer_save(request, cid):
	con = get_construct(request.user, cid)
	if request.method == 'POST' and con:
		data = json.loads(request.POST['data'])[0]
		con.pcrsettings.__dict__.update(data['pcrsettings'])
		con.pcrsettings.save()
		for f in data['fragments']:
			cf = ConstructFragment.objects.get(pk=f['id'])
			cf.concentration = f['t_c']
			cf.save()
			p = cf.primer_top()
			p.concentration = f['p_t_c']
			p.save()
			p = cf.primer_bottom()
			p.concentration = f['p_b_c']
			p.save()
		return HttpResponse('OK')
	else:
		return HttpResponseNotFound()