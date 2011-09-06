from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound
from pricalc.models import Oligo

from time import time
import json

def pricalc(request):
	t = loader.get_template('pricalc/pricalc.html')
	c = RequestContext(request, {
		'title':'Primer Calculator'
	})
	return HttpResponse(t.render(c))
	
def go(request):
	start = time()
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
			'SeqB':o.bottomPrimer(),
			'time':time()-start
		}),mimetype="application/json")