from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound

from digest.models import *


def calculator( request ):
	t = loader.get_template( 'tools/digest/calculator.html' )
	enzymes = Enzyme.objects.all()
	c = RequestContext( request, {
		'title':'Digest Calculator -> Calculator',
		'enzymes':enzymes,
	} )
	return HttpResponse( t.render( c ) )
	
def manufacturer_list( request ):
	t = loader.get_template( 'tools/digest/calculator_manufacturer_list.html' )
	enzyme_id = int( request.POST.get( 'enzyme' ) )
	if enzyme_id > 0:
		enzyme = Enzyme.objects.get( pk=enzyme_id )
		manufacturers = Manufacturer.objects.filter( brandedenzyme__in=enzyme.brandedenzyme_set.all() )
	else:
		manufacturers = [];
	c = RequestContext( request, {
		'manufacturers':manufacturers,
	} )
	return HttpResponse( t.render( c ) )
