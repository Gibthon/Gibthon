from django.db import models
from django.contrib.auth.models import User, UserManager
from annoying.fields import AutoOneToOneField
from django.contrib.contenttypes.models import ContentType

class GibthonUser(User):
	
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
		
	def constructmessage(self):
		return ConstructMessage.objects.filter(inbox=self)
	
	def fragmentmessage(self):
		return FragmentMessage.objects.filter(inbox=self)

class Message(models.Model):
	inbox = models.ForeignKey('Inbox', related_name='message')
	sender = models.CharField(max_length=50)
	data = models.TextField()
	received = models.DateField(auto_now_add=True)
	read = models.BooleanField()
	added = models.BooleanField()
	
	
class ConstructMessage(Message):

	def __unicode__(self):
		return "Message from " + self.sender + " to " + self.inbox.user.channel()

class FragmentMessage(Message):

	def __unicode__(self):
		return "Message from " + self.sender + " to " + self.inbox.user.channel()
	
