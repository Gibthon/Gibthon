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
	groups = BufferGroup.objects.all()
	c = RequestContext( request, {
		'title':'Digestion Calculator',
		'groups':groups,
	} )
	return HttpResponse( t.render( c ) )
	
def group( request ):
	t = loader.get_template( 'tools/digest/buffer_group.html' )
	groups = BufferGroup.objects.all()
	c = RequestContext( request, {
		'groups':groups,
	} )
	return HttpResponse( t.render( c ) )
	
def manufacturer( request ):
	t = loader.get_template( 'tools/digest/buffer_manufacturer.html' )
	ingredients = Ingredient.objects.all()
	buffers = [ RenderedBuffer( buffer, ingredients ) for buffer in Buffer.objects.all() ]
	c = RequestContext( request, {
		'ingredients':ingredients,
		'buffers':buffers,
	} )
	return HttpResponse( t.render( c ) )
	
def enzymes( request ):
	return HttpResponse()
