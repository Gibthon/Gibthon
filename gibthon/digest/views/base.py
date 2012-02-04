from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound

from digest.models import *

def base( request ):
	t = loader.get_template( 'tools/digest/base.html' )
	c = RequestContext( request, {
		'title':'Digestion Calculator',
	} )
	return HttpResponse( t.render( c ) )
	
	

