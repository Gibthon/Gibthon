from fragment.models import Gene, Feature, Qualifier, Annotation, Reference
from django.contrib import admin

admin.site.register(Qualifier)
admin.site.register(Gene)
admin.site.register(Feature)
admin.site.register(Annotation)
admin.site.register(Reference)