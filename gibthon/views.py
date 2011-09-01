from django.http import HttpResponse, HttpResponseRedirect
from django.template import Context, loader, RequestContext
from django.contrib.auth.decorators import login_required


import hashlib
import datetime

def redirect_home(request):
	return HttpResponseRedirect('/')

def molcal(request):
    return HttpResponseRedirect('/tools/molcalc/')

def ligcal(request):
    return HttpResponseRedirect('/tools/ligcalc/')
