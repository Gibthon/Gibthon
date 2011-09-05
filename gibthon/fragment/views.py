from fragment.models import *
from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound, Http404
from django.core.context_processors import csrf
from django.contrib.auth.decorators import login_required
from django.core.exceptions import *
from django.shortcuts import render_to_response

from Bio import SeqIO
from api import JsonResponse

genbank_mimetype = 'chemical/seq-na-genbank'
#genbank_mimetype = 'text/plain'

def get_fragment(user, fid):
	try:
		f = Gene.objects.get(pk=fid, owner=user)
	except ObjectDoesNotExist:
		return False
	else:
		return f

@login_required
def fragments(request):
	gene_list = Gene.objects.all().filter(owner=request.user)
	t = loader.get_template('fragment/fragments.html')
	c = RequestContext(request,{
		'fragment_list':gene_list,
		'title':'Fragments',
    })
	c.update(csrf(request))
	return HttpResponse(t.render(c))

@login_required	
def fragment(request, fid):
	f = get_fragment(request.user, fid)
	if f:
		f.origin = f.get_origin_display()
		t = loader.get_template('fragment/fragment.html')
		c = RequestContext(request,{
			'fragment':f,
			'title':'Fragment('+f.name+')',
		})
		return HttpResponse(t.render(c))
	else:
		t = loader.get_template('fragment/fragment_not_found.html')
		c = RequestContext(request,{
			'title':'Fragment not found',
		})
		return HttpResponse(t.render(c))	

@login_required
def fragment_meta(request, fid):
	t = loader.get_template('fragment/fragment_meta.html')
	c = RequestContext(request,{ 'id': fid,})
	return HttpResponse(t.render(c))

@login_required
def fragment_seq(request, fid):
	t = loader.get_template('fragment/sequence.html')
	c = RequestContext(request,{ 'id': fid,})
	return HttpResponse(t.render(c))

@login_required
def download(request, fid):
	f = get_fragment(request.user, fid)
	if f:
		response = HttpResponse(mimetype=genbank_mimetype)
		record = f.to_seq_record()
		SeqIO.write(record, response, 'genbank')
		return response
	else:
		raise Http404()

@login_required
def download_multi(request):
	"""download all the selected fragments in one file"""
	if request.method == 'POST' and 'selected' in request.POST:
		ids = request.POST.getlist('selected')
		genes = Gene.objects.filter(id__in = ids)
		records = [g.to_seq_record() for g in genes]
		response = HttpResponse(mimetype=genbank_mimetype)
		SeqIO.write(records, response, 'genbank')
		return response
	raise Http404()

@login_required
def delete(request):
	if request.method == 'POST':
		ids = []
		if 'selected' in request.POST:
			ids = request.POST.getlist('selected')
		elif 'selected[]' in request.POST:
			ids = request.POST.getlist('selected[]')
		else:
			raise Http404
			
		#remove the selected IDs from the database
		vids = []
		ok = []
		failed = []
		for id in ids:
			try:
				vids.append(int(id))
			except ValueError:
				failed.append(id)
		for id in vids:
			if Gene.remove(request.user, id):
				ok.append(id)
			else:
				failed.append(id)

		if not request.is_ajax():
			if len(failed) == 0:
				return HttpResponseRedirect('/fragment')
			return HttpResponseNotFound()
			
		return JsonResponse({'ok': ok, 'failed':failed,})

	raise Http404

