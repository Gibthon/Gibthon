from models import *
from forms import *
from fragment.models import Gene, Feature
from fragment.views import get_fragment

from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound, Http404
from django.core.context_processors import csrf
from django.contrib.auth.decorators import login_required
from django.views.decorators.http import condition
from django.core.exceptions import *
from django.conf import settings

from collections import OrderedDict
import csv, time, json, zipfile
from copy import copy
from cStringIO import StringIO
from xhtml2pdf import pisa


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
				'title': 'Construct Designer',
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
	if con and len(con.primer.all()) == 2*len(con.cf.all()) and len(con.primer.all()) > 0:
		t = loader.get_template('gibson/primers.html')
		c = RequestContext(request, {
			'construct': con,
			'primer_list': con.primer.all(),
			'title':'Construct Designer',
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
		cf_list = con.cf.all()
		feature_list = [FeatureListForm(cf, con) for cf in cf_list]
		list = zip(fragment_list, feature_list, cf_list)
		c = RequestContext(request,{
			'title':'Construct Designer',
			'list':list,
			'construct':con,
			'construct_form':construct_form,
		})
		return HttpResponse(t.render(c))
	else:
		return HttpResponseNotFound()

@login_required
def construct_fragment(request, cid):
	con = get_construct(request.user, cid)
	if con:
		t = loader.get_template('gibson/constructfragment.html')
		cf_list = con.cf.all()
		feature_list = [FeatureListForm(cf, con) for cf in cf_list]
		list = zip(cf_list, feature_list)
		c = RequestContext(request, {
			'list':list
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
		f = Gene.objects.all().filter(owner=request.user)
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
			con.add_fragment(f)
#			return HttpResponse('OK')
			return HttpResponseRedirect('/gibthon/'+cid+'/'+con.name+'/')
		else:
			return HttpResponseNotFound()
	else:
		return HttpResponseNotFound()

@login_required
def fragment_delete(request, cid, cfid):
	con = get_construct(request.user, cid)
	if con:
		try: cf = ConstructFragment.objects.get(id=cfid)
		except ObjectDoesNotExist: return HttpResponseNotFound()
		else:
			cf.delete()
			return HttpResponse("OK")
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

def pcr_step(name, temp, time):
	return {"type":"step", "name":name, "temp":str(temp), "time":str(time)}

def pcr_cycle(cf):
	tm = int((cf.primer_top().tm() + cf.primer_bottom().tm())/2 - 4)
	t = int(cf.fragment.length() * 45.0/1000.0)
	cycle = OrderedDict()
	cycle['type'] = 'cycle'
	cycle['count'] = 30
	cycle['steps'] = [pcr_step('Melting',98,10), pcr_step('Annealing', tm, 10), pcr_step('Elongation', 72, t)]
	pcr = OrderedDict()
	pcr['name'] = cf.construct.name + '-' + cf.fragment.name
	pcr['steps'] = [pcr_step('Initial Melting',98,30), cycle, pcr_step('Final Elongation', 72, 450),pcr_step('Final Hold', 4, 0)]
	pcr['lidtemp']  = 110
	return json.dumps(pcr, indent=4)

def pcr_instructions(request, cid):
	con = get_construct(request.user, cid)
	if con:
		response = HttpResponse(mimetype='application/zip')
		response['Content-Disposition'] = 'filename='+con.name+'-pcr.zip'
		pcr = [(con.name + '-' + cf.fragment.name+'.pcr',pcr_cycle(cf)) for cf in con.cf.all()]
		buffer = StringIO()
		zip = zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED)
		for name, f in pcr:
			zip.writestr(name, f)
		zip.close()
		buffer.flush()
		ret_zip = buffer.getvalue()
		buffer.close()
		response.write(ret_zip)
		return response
	else:
		return HttpReponseNotFound


@login_required
def primer_download(request, cid):
	con = get_construct(request.user, cid)
	if con:
		#set up response headers
		response = HttpResponse(mimetype='application/zip')
		response['Content-Disposition'] = 'attachment; filename='+con.name+'.zip'
		response.set_cookie('fileDownloadToken',request.GET['tk'])
		
		# get all the pcr instruction files
		pcr = [(con.name + '-' + cf.fragment.name+'.pcr',pcr_cycle(cf)) for cf in con.cf.all()]
		
		# write the csv file
		csvbuffer = StringIO()
		writer = csv.writer(csvbuffer)
		writer.writerow(['Name', 'Length', 'Melting Temperature', 'Sequence'])
		for p in con.primer.all():
			writer.writerow(p.csv())
		csvbuffer.flush()
		
		# write the pdf
		t = loader.get_template('gibson/pdf_primer.html')
		c = RequestContext(request,{
			'construct':con,
			'each':5.0/con.fragments.all().count()
		})
		pdfbuffer = StringIO()
		pdf = pisa.CreatePDF(StringIO(t.render(c).encode("ISO-8859-1")), pdfbuffer, link_callback=fetch_resources)
		
		# write the zip file
		zipbuffer = StringIO()
		zip = zipfile.ZipFile(zipbuffer, 'w', zipfile.ZIP_DEFLATED)
		# add the pcr files
		for name, f in pcr:
			zip.writestr(con.name+'/pcr/'+name, f)
		# add the csv file
		zip.writestr(con.name+'/primers.csv', csvbuffer.getvalue())
		# add the pdf
		zip.writestr(con.name+'/'+con.name+'.pdf', pdfbuffer.getvalue())
		# add the gb
		zip.writestr(con.name+'/'+con.name+'.gb', con.gb())
		# closing of buffers and return
		csvbuffer.close()
		pdfbuffer.close()
		zip.close()
		zipbuffer.flush()
		ret_zip = zipbuffer.getvalue()
		zipbuffer.close()
		response.write(ret_zip)
		return response
	else:
		return HttpResponseNotFound

@login_required
def pdf(request, cid):
	con = get_construct(request.user, cid)
	if con:
		response = HttpResponse(mimetype='application/pdf')
		t = loader.get_template('gibson/pdf_primer.html')
		c = RequestContext(request,{
			'construct':con,
			'each':5.0/con.fragments.all().count()
		})
		pdfbuffer = StringIO()
		pdf = pisa.CreatePDF(StringIO(t.render(c).encode("ISO-8859-1")), pdfbuffer, link_callback=fetch_resources)
		response.write(pdfbuffer.getvalue())
		pdfbuffer.close()
		return response
	else:
		return HttpResponseNotFound

def fetch_resources(uri, rel):
	path = os.path.join(settings.STATIC_ROOT, uri.replace(settings.STATIC_URL, "")[1:])
	return path
