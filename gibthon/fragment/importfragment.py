"""Contains JSON APIs to import fragments from various sources"""
from fragment.models import *
from django import forms
from django.contrib.auth.decorators import login_required
from django.conf.urls.defaults import patterns
from django.template import Context, loader, RequestContext
from django.core.context_processors import csrf
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound, Http404
from django.shortcuts import render_to_response

from api import JsonResponse, ERROR

@login_required
def fragment_import(request):
	"""return the import page"""
	t = loader.get_template('fragment/import.html')
	c = RequestContext(request,{})
	
	return HttpResponse(t.render(c))

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
	if request.method == 'GET' and 'id' in request.POST:
		id = request.GET['id']
		db = request.GET.get('database', 'nucleotide') #assume nucleotide by default
		handle = Entrez.efetch(db=db, id=id, rettype="gb")
		records = SeqIO.parse(handle, 'gb')
		for r in records:
			if len(Gene.objects.filter(name=r.name, owner=request.user)) == 0:
				Gene.add(r, 'NT', request.user)
		return JsonResponse("Imported id '%s' from Entrez." % id)
		
	raise Http404

entrezpatterns = patterns('',
	(r'^$', entrez),
	(r'^search/$', entrez_search),
	(r'^summary/$', entrez_summary),
	(r'^import/$', entrez_import),
)
