"""Contains JSON APIs to import fragments from various sources"""
from fragment.models import *
from django import forms
from django.contrib.auth.decorators import login_required
from django.conf.urls.defaults import patterns, include
from django.template import Context, loader, RequestContext
from django.core.context_processors import csrf
from django.core.files.uploadedfile import UploadedFile
from django.views.decorators.csrf import csrf_exempt
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound, Http404
from django.shortcuts import render_to_response

import simplejson as json
from gibthon.jsonresponses import JsonResponse, ERROR, RawJsonResponse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

@login_required
def fragment_import(request):
	"""return the import page"""
	t = loader.get_template('fragment/import.html')
	c = RequestContext(request,{})
	
	return HttpResponse(t.render(c))

##################################### MANUAL ############################

from forms import MetaForm, SequenceForm

@login_required
def manual_form(request):
	t = loader.get_template('fragment/MNform.html')
	meta = MetaForm()
	seq = SequenceForm()
	c = RequestContext(request, {	'meta': meta, 'seq': seq,
									'title':'Manually add a Fragment',
									'action': 'import/manual/add/',})
	
	return HttpResponse(t.render(c))

@login_required	
def manual_add(request):
	"""Add a fragment manually"""
	if request.method == 'GET':
		meta = MetaForm(request.GET)
		seq = SequenceForm(request.GET)
		if meta.is_valid() and seq.is_valid():
			record = SeqRecord(	seq.cleaned_data['seq'],
								name=meta.cleaned_data['name'], 
								id=meta.cleaned_data['name'],
								description=meta.cleaned_data['desc'])
			Gene.add(record, 'MN', request.user)
			return JsonResponse('OK');
		#figure out what the errors where
		errors = meta.errors
		errors.update(seq.errors)
		return JsonResponse(errors, ERROR)
	raise Http404

##################################### UPLOAD ############################


@login_required
def upload_form(request):
	"""return the import page"""
	t = loader.get_template('fragment/ULform.html')
	c = RequestContext(request,{})
	
	return HttpResponse(t.render(c))


fasta_types = ['fasta', 'fa', 'fsa', 'fna', 'ffn', 'faa', 'frn',]

@login_required
@csrf_exempt
def handle_upload(request):
	"""Handle a file upload"""
	if request.method == 'POST' and len(request.FILES) != 0:
		files = request.FILES.getlist('files[]')
		data = []
		for file in files:
			wrapped_file = UploadedFile(file)
			
			format = 'genbank' #assume genbank
			if wrapped_file.name.split('.')[-1] in fasta_types:
				format = 'fasta' #we were wrong!

			records = SeqIO.parse(wrapped_file, format)
			ids = []
			truncated = 0
			for record in records:
				(g, t) = Gene.add(record, 'UL', request.user, True)
				truncated += t['truncated']
				ids.append(g.id)
			
			file_dict = {	"name":wrapped_file.name, 
							"size":file.size,
							"error": None,
							"url":'/fragment/%s/' % ids[0],
							"delete_url":'/fragment/delete/', 
							"delete_type":"POST",
							"delete_data": json.dumps({'selected':ids,}),
						}
			
			if truncated > 0:
				e = ''
				if truncated == 1:
					e = "Warning: 1 annotation was truncated"
				else:
					e = "Warning: %i annotations were truncated" % truncated 
				file_dict['error'] = e
			
			data.append(file_dict)
		
		return RawJsonResponse(data)
	raise Http404
	
	
############################################################### PARTS #############################

from partsregistry import Part

@login_required
def part_form(request): #return the form
	t = loader.get_template('fragment/BBform.html')
	c = RequestContext(request, {})
	return HttpResponse(t.render(c))

@login_required
def part_import(request):
	"""API to import a part via AJAX"""
	if request.method == 'GET' and 'part' in request.GET:
		name = request.GET.get('part', '')
		try:
			part = Part(name)
			record = part.to_seq_record()
			Gene.add(record, 'BB', request.user)
		except Exception, e:
			return JsonResponse(e.message, ERROR)
		return JsonResponse('ok')
	raise Http404	
		

###################################################################################################
########## ENTREZ JSON API ########################################################################
###################################################################################################
########## Includes: search, summary and import ###################################################
###################################################################################################

entrez_databases = (	('nucleotide', 'All Nucleotide Databases'),
							('nuccore','Core Nucleotide'),
							('nucest','Expressed Sequence Tag (EST) Nucleotides'),
							('nucgss','Genome Survey Sequence (GSS) Nucleotides'),
							('gene','Gene'),
							('genome','Genome'),
						)

class EntrezSearchForm(forms.Form):
	database = forms.ChoiceField(choices=entrez_databases, label='Entrez Database')
	query = forms.CharField({'width': 60,}, label='Search Term')

@login_required
def entrez(request):
	"""return a blank entrez search form"""
	t = loader.get_template('fragment/NTform.html')
	c = RequestContext(request, {
			'fragment_form': EntrezSearchForm(),
			'type': 'NT',
			'action': 'search',})
	return HttpResponse(t.render(c))



@login_required
def entrez_search(request):
	"""	Handle JSON Entrez search request"""
	if request.method == 'GET':# and request.is_ajax():
		query = request.GET.get('query', '');
		db = request.GET.get('database', 'nucleotide') #assume nucleotide by default
		if query != '' and query != None:
			#fetch a list of Ids which match
			handle = Entrez.esearch(db=db, term=query)
			record = Entrez.read(handle)
			ids = record['IdList']
			if len(ids) == 0:
				return JsonResponse("No results for search '%s' on database '%s'." % (query, db), ERROR)
			else:
				return JsonResponse(ids)
		else:
			return JsonResponse("Error: Blank search query", ERROR);
		
	raise Http404

@login_required
def entrez_summary(request):
	"""	Handle JSON Entrez summary request"""
	if request.method == 'GET' and 'id' in request.GET:
		id = request.GET['id']
		db = request.GET.get('database', 'nucleotide') #assume nucleotide by default
		
		#fetch summary information for the id
		
		handle = Entrez.esummary(db=db, id=id)
		record = Entrez.read(handle)
		if record is None or len(record) < 1:
			return JsonResponse("Error: Could not get summary information for id '%s'" % id, ERROR)
		return JsonResponse(record[0]);	
	raise Http404

@login_required
def entrez_import(request):
	"""Handle a request to import an ID"""
	if request.method == 'GET' and 'id' in request.GET:
		id = request.GET['id']
		db = request.GET.get('database', 'nucleotide') #assume nucleotide by default
		handle = Entrez.efetch(db=db, id=id, rettype="gb")
		records = SeqIO.parse(handle, 'gb')
		for r in records:
			if len(Gene.objects.filter(name=r.name, owner=request.user)) == 0:
				Gene.add(r, 'NT', request.user)
		return JsonResponse("Imported id '%s' from Entrez." % id)
		
	raise Http404

