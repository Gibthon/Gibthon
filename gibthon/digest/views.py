from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound

from digest.models import *

def base( request ):
	t = loader.get_template( 'tools/digest/base.html' )
	c = RequestContext( request, {
		'title':'Digestion Calculator',
	} )
	return HttpResponse( t.render( c ) )
	
	
def buffers( request ):
	t = loader.get_template( 'tools/digest/buffers.html' )
	manufacturers = Manufacturer.objects.all()
	c = RequestContext( request, {
		'title':'Digestion Calculator',
		'manufacturers':manufacturers,
	} )
	return HttpResponse( t.render( c ) )
	
def group( request ):
	manufacturers = [ int(m) for m in request.POST.getlist('manufacturers[]') ]
	t = loader.get_template( 'tools/digest/buffer_group.html' )
	groups = [ buffer_group.get_filtered_buffers( manufacturers ) for buffer_group in BufferGroup.objects.all() ]
	c = RequestContext( request, {
		'groups':groups,
	} )
	return HttpResponse( t.render( c ) )
	
def manufacturer( request ):
	manufacturers = [ int(m) for m in request.POST.getlist('manufacturers[]') ]
	t = loader.get_template( 'tools/digest/buffer_manufacturer.html' )
	ingredients = Ingredient.objects.all()
	buffers = [ buffer.get_concentrations( ingredients ) for buffer in Buffer.objects.filter( manufacturer__id__in=manufacturers ) ]
	c = RequestContext( request, {
		'ingredients':ingredients,
		'buffers':buffers,
	} )
	return HttpResponse( t.render( c ) )
	
def enzymes( request ):
	return HttpResponse()
