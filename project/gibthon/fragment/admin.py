from fragment.models import Gene, Feature, Qualifier
from django.contrib import admin

admin.site.register(Qualifier)
admin.site.register(Gene)
admin.site.register(Feature)