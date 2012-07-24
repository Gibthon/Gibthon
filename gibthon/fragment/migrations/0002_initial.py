# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Gene'
        db.create_table('fragment_gene', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('owner', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=500)),
            ('sequence', self.gf('django.db.models.fields.TextField')(max_length=500000)),
            ('origin', self.gf('django.db.models.fields.CharField')(max_length=2)),
        ))
        db.send_create_signal('fragment', ['Gene'])

        # Adding model 'Reference'
        db.create_table('fragment_reference', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('gene', self.gf('django.db.models.fields.related.ForeignKey')(related_name='references', to=orm['fragment.Gene'])),
            ('title', self.gf('django.db.models.fields.CharField')(max_length=1024)),
            ('authors', self.gf('django.db.models.fields.CharField')(max_length=1024)),
            ('journal', self.gf('django.db.models.fields.CharField')(max_length=512)),
            ('medline_id', self.gf('django.db.models.fields.CharField')(max_length=24, blank=True)),
            ('pubmed_id', self.gf('django.db.models.fields.CharField')(max_length=24, blank=True)),
        ))
        db.send_create_signal('fragment', ['Reference'])

        # Adding model 'Annotation'
        db.create_table('fragment_annotation', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('gene', self.gf('django.db.models.fields.related.ForeignKey')(related_name='annotations', to=orm['fragment.Gene'])),
            ('key', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('value', self.gf('django.db.models.fields.CharField')(max_length=5120, blank=True)),
        ))
        db.send_create_signal('fragment', ['Annotation'])

        # Adding model 'Feature'
        db.create_table('fragment_feature', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('start', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('end', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('direction', self.gf('django.db.models.fields.CharField')(max_length=1)),
            ('gene', self.gf('django.db.models.fields.related.ForeignKey')(related_name='features', to=orm['fragment.Gene'])),
        ))
        db.send_create_signal('fragment', ['Feature'])

        # Adding model 'Qualifier'
        db.create_table('fragment_qualifier', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=30)),
            ('data', self.gf('django.db.models.fields.CharField')(max_length=512)),
            ('feature', self.gf('django.db.models.fields.related.ForeignKey')(related_name='qualifiers', to=orm['fragment.Feature'])),
        ))
        db.send_create_signal('fragment', ['Qualifier'])


    def backwards(self, orm):
        # Deleting model 'Gene'
        db.delete_table('fragment_gene')

        # Deleting model 'Reference'
        db.delete_table('fragment_reference')

        # Deleting model 'Annotation'
        db.delete_table('fragment_annotation')

        # Deleting model 'Feature'
        db.delete_table('fragment_feature')

        # Deleting model 'Qualifier'
        db.delete_table('fragment_qualifier')


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
        'fragment.annotation': {
            'Meta': {'object_name': 'Annotation'},
            'gene': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'annotations'", 'to': "orm['fragment.Gene']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'key': ('django.db.models.fields.CharField', [], {'max_length': '30'}),
            'value': ('django.db.models.fields.CharField', [], {'max_length': '5120', 'blank': 'True'})
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
        'fragment.qualifier': {
            'Meta': {'object_name': 'Qualifier'},
            'data': ('django.db.models.fields.CharField', [], {'max_length': '512'}),
            'feature': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'qualifiers'", 'to': "orm['fragment.Feature']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '30'})
        },
        'fragment.reference': {
            'Meta': {'object_name': 'Reference'},
            'authors': ('django.db.models.fields.CharField', [], {'max_length': '1024'}),
            'gene': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'references'", 'to': "orm['fragment.Gene']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'journal': ('django.db.models.fields.CharField', [], {'max_length': '512'}),
            'medline_id': ('django.db.models.fields.CharField', [], {'max_length': '24', 'blank': 'True'}),
            'pubmed_id': ('django.db.models.fields.CharField', [], {'max_length': '24', 'blank': 'True'}),
            'title': ('django.db.models.fields.CharField', [], {'max_length': '1024'})
        }
    }

    complete_apps = ['fragment']