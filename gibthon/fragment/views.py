from fragment.models import *
from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound, Http404
from django.core.context_processors import csrf
from django.contrib.auth.decorators import login_required
from django.core.exceptions import *
from django.shortcuts import render_to_response

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
		raise Http404()

@login_required
def add(request):
	
	if request.method == 'POST' and request.is_ajax():
		if request.POST['type'] == "BB":
			t = loader.get_template('fragment/BBform.html')
			form = BBForm()
			action = 'add'
		elif request.POST['type'] == "NT":
			t = loader.get_template('fragment/NTform.html')
			form = EntrezSearchForm()
			action = 'search'
		elif request.POST['type'] == "UL":
			t = loader.get_template('fragment/ULform.html')
			form = ULForm()
			action = 'add'
		c = RequestContext(request,{
			'fragment_form':form,
			'type': request.POST['type'],
			'action': action,
		})
		return HttpResponse(t.render(c))
	else:
		raise Http404()

entrez_databases = (	('nucleotide', 'All Nucleotide Databases'),
							('nuccore','Core Nucleotide'),
							('nucest','Expressed Sequence Tag (EST) Nucleotides'),
							('nucgss','Genome Survey Sequence (GSS) Nucleotides'),
							('gene','Gene'),
							('genome','Genome'),
						)

class EntrezSearchForm(forms.Form):
	database = forms.ChoiceField(choices=entrez_databases, label='Entrez Database')
	query = forms.CharField(label='Search Term')

@login_required
def entrez_search(request):
	"""Get search terms"""
	print "entrez_search(request): method='%s' is_ajax()='%s'" % (request.method, request.is_ajax())
	if request.method == 'POST' and request.is_ajax():
		print "entrez Search, post data: %s" % request.POST
		form = EntrezSearchForm(request.POST)
		if form.is_valid():
			#we have a valid search, time to make the search
			query = form.cleaned_data['query']
			database = form.cleaned_data['database']
			#fetch a list of Ids which match
			handle = Entrez.esearch(db=database, term=query)
			record = Entrez.read(handle)
			ids = record['IdList']
			#fetch summary information for each id
			summaries = []
			for id in ids:
				handle = Entrez.esummary(db=database, id=id)
				record = Entrez.read(handle)
				if record is None or len(record) < 1:
					continue
				summaries.append(record[0])

			return render_to_response("fragment/result_form.html",{
												'action': 'add',
												'type': 'NT',
												'summaries': summaries,
												'database': database,},
												context_instance=RequestContext(request))
	raise Http404()


@login_required
def search(request, type):
	if request.method == 'POST':
		if type == "BB":
			return HttpResponseNotFound()
		elif type == "NT":
			print 'got search type NT'
			return entrez_search(request)
		elif type == "UL":
			return HttpResponseNotFound()
		return HttpResponseRedirect('/fragment')
	return HttpResponseNotFound()

@login_required
def entrez_import(request):
	if request.method == 'POST' and 'import' in request.POST and 'database' in request.POST:
		ids = request.POST.getlist('import')
		#TODO, check that we were handed back valid IDs
		database = request.POST['database']
		for id in ids:
			handle = Entrez.efetch(db=database, id=id, rettype="gb")
			records = SeqIO.parse(handle, 'gb')
			for r in records:
				if len(Gene.objects.filter(name=r.name, owner=request.user)) == 0:
					Gene.add(r, 'NT', request.user)


@login_required
def add_submit(request, type):
	if request.method == 'POST':
		if type == "BB":
			form = BBForm(request.POST)
			if form.is_valid():
				#if len(Gene.objects.filter(name=form.cleaned_data['id'], owner=request.user)) == 0:
				#	Gene.add(form.cleaned_data['id'], type, request.user)
				print 'partsregistry import not implemented!'
		elif type == "NT":
			entrez_import(request)
			return HttpResponseRedirect('/fragment')
		elif type == "UL":
			form = ULForm(request.POST, request.FILES)
			if form.is_valid():
				records = SeqIO.parse(form.cleaned_data['file'])
				for record in records:
					if len(Gene.objects.filter(name=record.name, owner=request.user)) == 0:
						Gene.add(record, type, request.user)
		return HttpResponseRedirect('/fragment')
	return HttpResponseNotFound()