from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound

from digest.models import *


	
def enzymes( request ):
	t = loader.get_template( 'tools/digest/enzymes.html' )
	manufacturers = Manufacturer.objects.all()
	c = RequestContext( request, {
		'title':'Digestion Calculator -> Enzymes',
		'manufacturers':manufacturers,
	} )
	return HttpResponse( t.render( c ) )
	
def manufacturer( request ):
	t = loader.get_template( 'tools/digest/enzymes_manufacturer.html' )
	manufacturer = Manufacturer.objects.get( pk=int( request.POST.get('manufacturer' ) ) )
	enzymes = Enzyme.objects.filter( manufacturer=manufacturer )
	c = RequestContext( request, {
		'manufacturers':manufacturers,
		'enzymes':enzymes,
	} )
	return HttpResponse( t.render( c ) )
