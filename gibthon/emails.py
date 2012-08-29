from django.core.mail import EmailMessage
from django.conf import settings
import hashlib
import datetime


def RegisterEmail1(to):
	email_hash = hashlib.md5(settings.SECRET_KEY + to + str(datetime.date.today())).hexdigest()
	email = EmailMessage('Registration at Gibthon.org',
'''Thank you for registering at gibthon.org.
				
Please click on the following link to create your account.
				
http://www.gibthon.org/user/register/%s/

Using the email %s

This link will expire at midnight.

Regards,

Gibthon Admin'''%(email_hash, to),
		'admin@gibthon.org',
		[to],
		['users@gibthon.org']
	)
	email.send()


def RegisterEmail2(user):
	email = EmailMessage('Registration at Gibthon.org',
'''Dear %s %s,

Thank you for registering at gibthon.org.
				
Your account has now been activated

Regards,

Gibthon Admin'''%(user.first_name, user.last_name),
		'admin@gibthon.org',
		[user.email],
		['users@gibthon.org']
	)
	email.send()
