from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound

from digest.models import *

def digest( request ):
	t = loader.get_template( 'tools/digest/base.html' )
	ingredients = Ingredient.objects.all()
	manufacturers = Manufacturer.objects.all()
	c = RequestContext( request, {
		'title':'Digestion Calculator',
		'ingredients':ingredients,
		'manufacturers':manufacturers,
	} )
	return HttpResponse( t.render( c ) )
