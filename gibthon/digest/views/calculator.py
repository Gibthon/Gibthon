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
	enzyme = BrandedEnzyme.objects.get( pk=request.POST.get('enzyme') )
	buffers = enzyme.enzymereactionbuffer_set.all()
	c = RequestContext( request, {
		'enzyme':enzyme,
		'buffers':buffers,
	} )
	return HttpResponse( t.render( c ) )
	
def get_buffer_double( request ):
	brandedenzymes = [ BrandedEnzyme.objects.get( pk=eid ) for eid in request.POST.getlist('enzymes[]') ]

	if brandedenzymes[0].manufacturer == brandedenzymes[1].manufacturer:	
		t = loader.get_template( 'tools/digest/calculator_double_enzyme.html' )
		buffers = brandedenzymes[0].buffers.all()
		print buffers
		c = RequestContext( request, {
			'brandedenzymes':brandedenzymes,
			'buffers':buffers,
		} )
	else:
		t = loader.get_template( 'tools/digest/calculator_double_trouble_enzyme.html' )
		enzyme1 = brandedenzymes[0]
		enzyme2 = brandedenzymes[1]
		enzyme2in1 = enzyme2.manufacturer.brandedenzyme_set.filter( enzyme=enzyme1.enzyme ).get()
		enzyme1in2 = enzyme1.manufacturer.brandedenzyme_set.filter( enzyme=enzyme2.enzyme ).get()
		buffers1 = enzyme1.buffers.all()
		buffers2 = enzyme2.buffers.all()
		enzyme2.get_related_activities( buffers1 )
		enzyme2.get_activities( buffers2 )
		enzyme1.get_related_activities( buffers2 )
		enzyme1.get_activities( buffers1 )
		enzyme2in1.get_activities( buffers2 )
		enzyme1in2.get_activities( buffers1 )
		c = RequestContext( request, 
			{ 'enzyme1':enzyme1,
			  'enzyme2':enzyme2,
			  'buffers1':buffers1,
			  'buffers2':buffers2,
			  'enzyme2in1':enzyme2in1,
			  'enzyme1in2':enzyme1in2,
			} )
	return HttpResponse( t.render( c ) )
