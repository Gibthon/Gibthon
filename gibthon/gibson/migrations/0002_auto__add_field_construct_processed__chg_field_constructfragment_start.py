# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding field 'Construct.processed'
        db.add_column('gibson_construct', 'processed',
                      self.gf('django.db.models.fields.BooleanField')(default=False),
                      keep_default=False)


        # Changing field 'ConstructFragment.start_feature'
        db.alter_column('gibson_constructfragment', 'start_feature_id', self.gf('django.db.models.fields.related.ForeignKey')(null=True, to=orm['fragment.Feature']))

        # Changing field 'ConstructFragment.end_feature'
        db.alter_column('gibson_constructfragment', 'end_feature_id', self.gf('django.db.models.fields.related.ForeignKey')(null=True, to=orm['fragment.Feature']))

    def backwards(self, orm):
        # Deleting field 'Construct.processed'
        db.delete_column('gibson_construct', 'processed')


        # User chose to not deal with backwards NULL issues for 'ConstructFragment.start_feature'
        raise RuntimeError("Cannot reverse this migration. 'ConstructFragment.start_feature' and its values cannot be restored.")

        # User chose to not deal with backwards NULL issues for 'ConstructFragment.end_feature'
        raise RuntimeError("Cannot reverse this migration. 'ConstructFragment.end_feature' and its values cannot be restored.")

    models = {
        'auth.group': {
            'Meta': {'object_name': 'Group'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        'auth.permission': {
            'Meta': {'ordering': "('content_type__app_label', 'content_type__model', 'codename')", 'unique_together': "(('content_type', 'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['contenttypes.ContentType']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.Group']", 'symmetrical': 'False', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        'fragment.feature': {
            'Meta': {'ordering': "['start']", 'object_name': 'Feature'},
            'direction': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'end': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'gene': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'features'", 'to': "orm['fragment.Gene']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'start': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '30'})
        },
        'fragment.gene': {
            'Meta': {'object_name': 'Gene'},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '500'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'origin': ('django.db.models.fields.CharField', [], {'max_length': '2'}),
            'owner': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['auth.User']"}),
            'sequence': ('django.db.models.fields.TextField', [], {'max_length': '500000'})
        },
        'gibson.construct': {
            'Meta': {'object_name': 'Construct'},
            'created': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.CharField', [], {'max_length': '2000'}),
            'fragments': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "'construct_slave'", 'blank': 'True', 'through': "orm['gibson.ConstructFragment']", 'to': "orm['fragment.Gene']"}),
            'genbank': ('django.db.models.fields.related.OneToOneField', [], {'blank': 'True', 'related_name': "'construct_master'", 'unique': 'True', 'null': 'True', 'to': "orm['fragment.Gene']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modified': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '80'}),
            'owner': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['auth.User']", 'null': 'True'}),
            'processed': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'shape': ('django.db.models.fields.CharField', [], {'max_length': '1'})
        },
        'gibson.constructfragment': {
            'Meta': {'ordering': "['order']", 'object_name': 'ConstructFragment'},
            'concentration': ('django.db.models.fields.DecimalField', [], {'default': '100', 'max_digits': '4', 'decimal_places': '1'}),
            'construct': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'cf'", 'to': "orm['gibson.Construct']"}),
            'direction': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'end_feature': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'end_feature'", 'null': 'True', 'to': "orm['fragment.Feature']"}),
            'end_offset': ('django.db.models.fields.IntegerField', [], {}),
            'fragment': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'cf'", 'to': "orm['fragment.Gene']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'order': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'start_feature': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'start_feature'", 'null': 'True', 'to': "orm['fragment.Feature']"}),
            'start_offset': ('django.db.models.fields.IntegerField', [], {})
        },
        'gibson.pcrsettings': {
            'Meta': {'object_name': 'PCRSettings'},
            'buffer_d': ('django.db.models.fields.DecimalField', [], {'default': '1', 'max_digits': '3', 'decimal_places': '1'}),
            'buffer_s': ('django.db.models.fields.DecimalField', [], {'default': '10', 'max_digits': '3', 'decimal_places': '1'}),
            'construct': ('annoying.fields.AutoOneToOneField', [], {'related_name': "'pcrsettings'", 'unique': 'True', 'to': "orm['gibson.Construct']"}),
            'dntp_d': ('django.db.models.fields.DecimalField', [], {'default': '0.8', 'max_digits': '3', 'decimal_places': '1'}),
            'dntp_s': ('django.db.models.fields.DecimalField', [], {'default': '10', 'max_digits': '3', 'decimal_places': '1'}),
            'enzyme_d': ('django.db.models.fields.DecimalField', [], {'default': '2.5', 'max_digits': '3', 'decimal_places': '1'}),
            'enzyme_s': ('django.db.models.fields.DecimalField', [], {'default': '2.5', 'max_digits': '3', 'decimal_places': '1'}),
            'error_margin': ('django.db.models.fields.PositiveSmallIntegerField', [], {'default': '10'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'primer_d': ('django.db.models.fields.DecimalField', [], {'default': '0.4', 'max_digits': '3', 'decimal_places': '1'}),
            'repeats': ('django.db.models.fields.PositiveSmallIntegerField', [], {'default': '1'}),
            'template_d': ('django.db.models.fields.DecimalField', [], {'default': '100', 'max_digits': '4', 'decimal_places': '1'}),
            'volume_each': ('django.db.models.fields.DecimalField', [], {'default': '12.5', 'max_digits': '3', 'decimal_places': '1'})
        },
        'gibson.primer': {
            'Meta': {'ordering': "['stick']", 'object_name': 'Primer'},
            'boxplot': ('django.db.models.fields.files.ImageField', [], {'max_length': '100'}),
            'concentration': ('django.db.models.fields.DecimalField', [], {'default': '5', 'max_digits': '4', 'decimal_places': '1'}),
            'construct': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'primer'", 'to': "orm['gibson.Construct']"}),
            'flap': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "'flap'", 'unique': 'True', 'to': "orm['gibson.PrimerHalf']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '80'}),
            'stick': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "'stick'", 'unique': 'True', 'to': "orm['gibson.PrimerHalf']"})
        },
        'gibson.primerhalf': {
            'Meta': {'object_name': 'PrimerHalf'},
            'cfragment': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'ph'", 'to': "orm['gibson.ConstructFragment']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'length': ('django.db.models.fields.PositiveSmallIntegerField', [], {}),
            'top': ('django.db.models.fields.BooleanField', [], {'default': 'False'})
        },
        'gibson.settings': {
            'Meta': {'object_name': 'Settings'},
            'construct': ('annoying.fields.AutoOneToOneField', [], {'related_name': "'settings'", 'unique': 'True', 'to': "orm['gibson.Construct']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'mg_salt': ('django.db.models.fields.DecimalField', [], {'default': '0.0', 'max_digits': '3', 'decimal_places': '2'}),
            'min_anneal_tm': ('django.db.models.fields.PositiveSmallIntegerField', [], {'default': '50'}),
            'min_overlap': ('django.db.models.fields.PositiveSmallIntegerField', [], {'default': '20'}),
            'min_primer_tm': ('django.db.models.fields.PositiveSmallIntegerField', [], {'default': '60'}),
            'na_salt': ('django.db.models.fields.DecimalField', [], {'default': '1', 'max_digits': '3', 'decimal_places': '2'}),
            'ss_safety': ('django.db.models.fields.PositiveSmallIntegerField', [], {'default': '3'})
        },
        'gibson.warning': {
            'Meta': {'object_name': 'Warning'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'primer': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'warning'", 'to': "orm['gibson.Primer']"}),
            'text': ('django.db.models.fields.CharField', [], {'max_length': '150'}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '2'})
        }
    }

    complete_apps = ['gibson']