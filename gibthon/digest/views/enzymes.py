from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound

from digest.models import *


	
def enzymes( request ):
	t = loader.get_template( 'tools/digest/enzymes.html' )
	c = RequestContext( request, {
		'title':'Digestion Calculator -> Enzymes',
		'manufacturers':Manufacturer.objects.all(),
		'enzymes':Enzyme.objects.all(),
	} )
	return HttpResponse( t.render( c ) )
	
def manufacturer( request ):
	t = loader.get_template( 'tools/digest/enzymes_manufacturer.html' )
	manufacturer = Manufacturer.objects.get( pk=int( request.POST.get('manufacturer' ) ) )
	print request.POST.getlist('enzymes[]')
	enzymes = BrandedEnzyme.objects.filter( manufacturer=manufacturer).filter( enzyme__id__in=request.POST.getlist('enzymes[]') )
	c = RequestContext( request, {
		'manufacturer':manufacturer,
		'enzymes':enzymes,
	} )
	return HttpResponse( t.render( c ) )
