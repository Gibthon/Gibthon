from molcal.models import *

from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound
from django.core.context_processors import csrf
from django.contrib.auth.decorators import login_required
from django.views.decorators.http import condition
from django.core.exceptions import *

import simplejson as json

def molcalc(request):
	t = loader.get_template("tools/molcal.html")
	c = RequestContext(request, {
		'title':'Molarity Calculator',
	})
	return HttpResponse(t.render(c))
	
def getall(request):
	term = request.GET['term']
	mols = Mol.objects.filter(name__istartswith=term)
	json_data = json.dumps([{'value':m.name, 'label':m.name, 'fullname':m.fullname, 'formula':m.formula, 'weight':str(m.weight)} for m in mols])
	return HttpResponse(json_data)