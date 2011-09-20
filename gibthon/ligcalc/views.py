from fragment.models import Gene
import simplejson as json
from django.http import HttpResponse
from django.contrib.auth.decorators import login_required
from django.template import Context, loader, RequestContext


def ligcalc(request):
	t = loader.get_template("tools/ligate.html")
	c = RequestContext(request, {
		'title':'Ligation Calculator',
	})
	return HttpResponse(t.render(c))

@login_required
def get(request):
	term = request.GET['term']
	fragments = Gene.objects.filter(name__icontains=term, owner=request.user)
	json_data = json.dumps([{'value':f.name, 'label':f.name, 'length':f.length()} for f in fragments])
	return HttpResponse(json_data)