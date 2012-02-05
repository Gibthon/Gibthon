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
		brandedenzymes = BrandedEnzyme.objects.filter( enzyme=enzyme )
	else:
		brandedenzymes = [];
	c = RequestContext( request, {
		'brandedenzymes':brandedenzymes,
	} )
	return HttpResponse( t.render( c ) )
	
def get_buffer_single( request ):
	t = loader.get_template( 'tools/digest/calculator_single_enzyme.html' )
	brandedenzyme = BrandedEnzyme.objects.get( pk=request.POST.get('enzyme') )
	enzymereactionbuffers = EnzymeReactionBuffer.objects.filter( enzyme=brandedenzyme )
	c = RequestContext( request, {
		'brandedenzyme':brandedenzyme,
		'enzymereactionbuffers':enzymereactionbuffers,
	} )
	return HttpResponse( t.render( c ) )
	
def get_buffer_double( request ):
	t = loader.get_template( 'tools/digest/calculator_double_enzyme.html' )
	brandedenzymes = [ BrandedEnzyme.objects.get( pk=eid ) for eid in request.POST.getlist('enzymes[]') ]
	buffers = brandedenzymes[0].buffers.all()
	c = RequestContext( request, {
		'brandedenzymes':brandedenzymes,
		'buffers':buffers,
	} )
	return HttpResponse( t.render( c ) )
