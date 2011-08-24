from django.db import models
from django.contrib.auth.models import User, UserManager
from annoying.fields import AutoOneToOneField
from django.contrib.contenttypes.models import ContentType

from gibthon.messages import MessagePasser

import json

class GibthonUser(User):

	channel_key = models.CharField(max_length=120)
	objects = UserManager()
	
	def channel(self):
		return 'GCD-' + self.username

class Inbox(models.Model):
	user = AutoOneToOneField('GibthonUser', related_name='inbox')
	
	def __unicode__(self):
		return self.user.username + "'s inbox"
	
	def unread(self):
		return self.message.filter(read=False)
		
	def not_added(self):
		return self.message.filter(added=False)
		
	def messagePasser(self):
		return MessagePasser(self.user.channel())
	
	def fetch(self):
		messages = self.messagePasser().fetch()
		if messages == []:
			return True
		else:
			for message in messages:
				if message == '':
					continue
				_data = json.dumps(message['data'])
				_sender = message['source']
				_origin = message['source'].lower()
				try:
					message['type']
				except:
					_type = 'cn'
				else:
					_type = message['type']
				m = Message.objects.create(inbox=self, sender=_sender, data=_data, type=_type, origin=_origin)
			self.messagePasser().clear()
				

class Message(models.Model):
	inbox = models.ForeignKey('Inbox', related_name='message')
	sender = models.CharField(max_length=50)
	data = models.TextField()
	received = models.DateTimeField(auto_now_add=True)
	read = models.BooleanField(default=False)
	added = models.BooleanField(default=False)
	TYPE_CHOICES = (
		('cn', 'Construct'),
		('fr', 'Fragment'),
	)
	type = models.CharField(max_length=2, choices=TYPE_CHOICES)
	ORIGIN_CHOICES = (
		('pr', 'Parts Registry'),
		('gcd', 'Gibthon Construct Designer'),
		('gec', 'Genetic Engineering of Cells'),
		('xxx', 'Unknown'),
	)
	origin = models.CharField(max_length=3, choices=ORIGIN_CHOICES)
	
	def description(self):
		return self.get_type_display() + ' from ' + self.get_origin_display()
	
	def __unicode__(self):
		return "Message from " + self.sender + " to " + self.inbox.user.channel()
