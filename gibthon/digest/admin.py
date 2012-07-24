from digest.models import Manufacturer, Enzyme, Buffer, Ingredient, BufferIngredient, BrandedEnzyme, EnzymeReactionBuffer, BufferGroup
from django.contrib import admin

admin.site.register( Manufacturer )
admin.site.register( Enzyme )
admin.site.register( Buffer )
admin.site.register( Ingredient )
admin.site.register( BufferIngredient )
admin.site.register( BrandedEnzyme )
admin.site.register( EnzymeReactionBuffer )
admin.site.register( BufferGroup )
