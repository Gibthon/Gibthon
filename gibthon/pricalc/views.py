from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound

def pricalc(request):
	t = loader.get_template('pricalc/pricalc.html')
	c = RequestContext(request, {
		'title':'Primer Calculator'
	})
	return HttpResponse(t.render(c))
	
def go(request):
	print request.POST
	left = request.POST.get('gene1')
	right = request.POST.get('gene2')
	lefthalf = left[::-1][int(request.POST.get('sl'))-1:int(request.POST.get('el'))][::-1]
	righthalf = right[int(request.POST.get('sr'))-1:int(request.POST.get('er'))]
	primer = lefthalf+righthalf
	return HttpResponse(primer[::-1])