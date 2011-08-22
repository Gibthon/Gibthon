from fragment.models import *
from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound
from django.core.context_processors import csrf
from django.contrib.auth.decorators import login_required
from django.core.exceptions import *

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
def download(request, fid):
	f = get_fragment(request.user, fid)
	if f:
		#response = HttpResponse(mimetype='chemical/seq-na-genbank')
		response = HttpResponse(mimetype='text/plain')
		response.write(f.gb())
		return response
	else:
		return HttpResponseNotFound()

@login_required
def add(request):
	if request.method == 'POST' and request.is_ajax():
		if request.POST['type'] == "BB":
			t = loader.get_template('fragment/BBform.html')
			form = BBForm()
		elif request.POST['type'] == "NT":
			t = loader.get_template('fragment/NTform.html')
			form = NTForm()
		elif request.POST['type'] == "UL":
			t = loader.get_template('fragment/ULform.html')
			form = ULForm()
		c = RequestContext(request,{
			'fragment_form':form,
			'type':request.POST['type'],
		})
		return HttpResponse(t.render(c))
	else:
		return HttpResponseNotFound()

@login_required
def add_submit(request, type):
	if request.method == 'POST':
		ok = False
		if type == "BB":
			form = BBForm(request.POST)
			if form.is_valid():
				if len(Gene.objects.filter(name=form.cleaned_data['id'], owner=request.user)) == 0:
					ok = Gene.add(form.cleaned_data['id'], type, request.user)
		elif type == "NT":
			form = NTForm(request.POST)
			if form.is_valid():
				if len(Gene.objects.filter(name=form.cleaned_data['id'], owner=request.user)) == 0:
					ok = Gene.add(form.cleaned_data['id'], type, request.user)
		elif type == "UL":
			form = ULForm(request.POST, request.FILES)
			if form.is_valid():
				if len(Gene.objects.filter(name=form.cleaned_data['file'].name.split('.')[0], owner=request.user)) == 0:
					ok = Gene.add(form.cleaned_data['file'], type, request.user)
		if ok:
			return HttpResponseRedirect('/fragment')
		else:
			return HttpResponseNotFound()
	return HttpResponseNotFound()