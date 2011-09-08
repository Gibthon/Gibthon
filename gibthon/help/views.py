# Create your views here.

from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound

def help(request):
	return HttpResponseRedirect('/')

def gibson(request):
	t = loader.get_template('help/gibson.html')
	c = RequestContext(request, {
		'title':'Help -> Gibson'
	})
	return HttpResponse(t.render(c))
	
def primer(request):
	t = loader.get_template('help/primer.html')
	c = RequestContext(request, {
		'title':'Help -> Primer Design'
	})
	return HttpResponse(t.render(c))

def gibthon(request):
	t = loader.get_template('help/gibthon.html')
	c = RequestContext(request, {
		'title':'Help -> Construct Designer'
	})
	return HttpResponse(t.render(c))