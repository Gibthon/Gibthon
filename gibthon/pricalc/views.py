from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound
from pricalc.models import Oligo
from django.contrib.auth.decorators import login_required


from fragment.models import Gene

import json

def pricalc(request):
	t = loader.get_template('pricalc/pricalc.html')
	c = RequestContext(request, {
		'title':'Primer Calculator'
	})
	return HttpResponse(t.render(c))
	
def go(request):
	try:
		o = Oligo(data = request.POST)
	except Exception as e:
		print e
	return HttpResponse(
		json.dumps({
			'TmT':o.topTm(),
			'TmB':o.bottomTm(),
			'TmF':o.fullTm(),
			'SeqT':o.topPrimer(),
			'SeqB':o.bottomPrimer()
		}),mimetype="application/json")
		
def selfprime(request):
	try:
		o = Oligo(data = request.POST)
	except Exception as e:
		print e
	image, warnings = o.selfPrime()
	return HttpResponse(
		json.dumps({
			'warnings':warnings,
			'image':image
		}), mimetype="application/json")
		
@login_required
def list(request):
	term = request.GET['term']
	fragments = Gene.objects.filter(name__icontains=term, owner=request.user)
	json_data = json.dumps([{'value':f.name, 'label':f.name, 'first':f.sequence[:50], 'last':f.sequence[-50:]} for f in fragments])
	return HttpResponse(json_data)