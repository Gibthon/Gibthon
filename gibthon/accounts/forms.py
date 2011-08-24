from django import forms
from django.contrib.auth.models import User
from django.contrib.auth import authenticate
from gibthon.accounts.models import GibthonUser
from captcha.fields import CaptchaField
from django.core.mail import EmailMessage
from django.conf import settings
import hashlib, datetime


class UserRegisterForm1(forms.Form):
	email = forms.EmailField(label=("E-mail"), max_length=75)
	captcha = CaptchaField()
	
	def clean_email(self):
		email = self.cleaned_data["email"]
		try:
			User.objects.get(email=email)
		except User.DoesNotExist:
			return email
		raise forms.ValidationError("Someone is already using that email address")

class UserRegisterForm2(forms.ModelForm):
	"""
	A form that creates a user, with no privileges, from the given username and password.
	"""
	username = forms.RegexField(label=("Username"), max_length=30, regex=r'^[\w.@+-]+$',
		help_text = ("Required. 30 characters or fewer. Letters, digits and @/./+/-/_ only."),
		error_messages = {'invalid': ("This value may contain only letters, numbers and @/./+/-/_ characters.")},)
	password1 = forms.CharField(label=("Password"), widget=forms.PasswordInput)
	password2 = forms.CharField(label=("Password confirmation"), widget=forms.PasswordInput,
		help_text = ("Enter the same password as above, for verification."))

	class Meta:
		model = GibthonUser
		fields = ("username","email","first_name", "last_name")

	def is_valid(self, hash, *args, **kwargs):
		self.hash = hash
		return super(UserRegisterForm2, self).is_valid(*args, **kwargs)

	def clean_first_name(self):
		first_name = self.cleaned_data["first_name"]
		if first_name != '':
			return first_name
		raise forms.ValidationError(("Please enter your first name."))

	def clean_last_name(self):
		last_name = self.cleaned_data["last_name"]
		if last_name != '':
			return last_name
		raise forms.ValidationError(("Please enter your last name."))

	def clean_username(self):
		username = self.cleaned_data["username"]
		try:
			User.objects.get(username=username)
		except User.DoesNotExist:
			return username
		raise forms.ValidationError(("A user with that username already exists."))

	def clean_email(self):
		email = self.cleaned_data["email"]
		try:
			User.objects.get(email=email)
		except User.DoesNotExist:
			if self.hash != hashlib.md5(settings.SECRET_KEY + email + str(datetime.date.today())).hexdigest():
				raise forms.ValidationError(("Please use the email address you used to get this URL"))
			return email
		raise forms.ValidationError(("That email has already been used."))

	def clean_password2(self):
		password1 = self.cleaned_data.get("password1", "")
		password2 = self.cleaned_data["password2"]
		if password1 != password2:
			raise forms.ValidationError(("The two password fields didn't match."))
		return password2

	def save(self, commit=True):
		user = super(UserRegisterForm2, self).save(commit=False)
		user.set_password(self.cleaned_data["password1"])
		if commit:
			user.save()
		user = authenticate(username = user.username, password=self.cleaned_data["password1"])
		email = EmailMessage('Registration at Gibthon.org',
'''Thank you for registering at gibthon.org.
				
Your account has now been activated

Regards,

Gibthon Admin''',
			'admin@gibthon.org',
			[self.cleaned_data['email']],
			['users@gibthon.org'],)
		email.send()
		return user
